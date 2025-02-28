import pickle
import atexit
import sys
import uuid
from multiprocessing import shared_memory
from typing import Dict, Any, List, Optional


# Global dictionary to track all shared memory blocks
_shared_memory_blocks = {}


# Function to clean up shared memory resources on program exit
def cleanup_shared_memory():
    """Clean up all shared memory objects when the program exits."""
    global _shared_memory_blocks
    
    if _shared_memory_blocks:
        try:
            # Close and unlink all shared memory blocks
            for name, shm in list(_shared_memory_blocks.items()):
                try:
                    shm.close()
                    shm.unlink()
                    # No need to print anything during normal exit
                except Exception as e:
                    if not sys.stdout.closed:
                        print(f"Error cleaning up shared memory '{name}': {e}")
            
            # Clear the dictionary
            _shared_memory_blocks.clear()
        except Exception as e:
            if not sys.stdout.closed:
                print(f"Error during shared memory cleanup: {e}")


# Register the cleanup function to run on program exit
atexit.register(cleanup_shared_memory)


class SharedDict:
    """A dictionary-like object that provides fast read access to shared memory data.
    
    This class provides dictionary-like access to shared memory with performance
    approaching that of a local dictionary. On first access in each worker process,
    the dictionary is loaded into local memory for maximum lookup performance while
    still saving memory across the main process.
    """
    
    def __init__(self, shared_dict_info, dict_id=None):
        """Initialize with shared memory information.
        
        Args:
            shared_dict_info (dict): Information about the shared memory block
            dict_id (str, optional): Unique identifier for this shared dictionary
        """
        self.info = shared_dict_info
        self.dict_id = dict_id
        # Initialize with an empty dictionary to satisfy type checkers
        self._dict = {}
        self._shm = None
        self._loaded = False
    
    def _load_dict(self):
        """Load the dictionary from shared memory if not already loaded.
        
        This only happens once per process, making all subsequent lookups
        as fast as a local dictionary.
        """
        if not self._loaded:
            try:
                # Access the shared memory block
                shm = shared_memory.SharedMemory(name=self.info['name'])
                # Deserialize the dictionary
                data = bytes(shm.buf[:self.info['size']])
                self._dict = pickle.loads(data)
                # Keep a reference to the shared memory to prevent it from being garbage collected
                self._shm = shm
                self._loaded = True
            except Exception as e:
                print(f"Error loading shared dictionary: {e}")
                # Ensure we have a valid dictionary even if loading fails
                if not self._dict:
                    self._dict = {}
    
    def __del__(self):
        """Clean up shared memory when this object is garbage collected.
        
        This ensures that shared memory is released in worker processes.
        Note: The main process still needs the atexit handler to clean up the
        original shared memory block.
        """
        if hasattr(self, '_shm') and self._shm is not None:
            try:
                # Only close it - don't unlink as the main process manages that
                self._shm.close()
                self._shm = None
            except Exception:
                # Silently ignore errors during garbage collection
                pass
    
    def get(self, key, default=None):
        """Get a value for key, or default if key is not present.
        
        Args:
            key: The key to look up
            default: Value to return if key is not found
            
        Returns:
            The value for key if present, otherwise default
        """
        self._load_dict()
        return self._dict.get(key, default)
    
    def __getitem__(self, key):
        """Get a value for key.
        
        Args:
            key: The key to look up
            
        Returns:
            The value for key
            
        Raises:
            KeyError: If key is not found
        """
        self._load_dict()
        return self._dict[key]
    
    def __contains__(self, key):
        """Check if key is in the dictionary.
        
        Args:
            key: The key to check for
            
        Returns:
            bool: True if key is present, False otherwise
        """
        self._load_dict()
        return key in self._dict
    
    def keys(self):
        """Get the dictionary keys.
        
        Returns:
            An iterable of keys
        """
        self._load_dict()
        return self._dict.keys()
    
    def values(self):
        """Get the dictionary values.
        
        Returns:
            An iterable of values
        """
        self._load_dict()
        return self._dict.values()
    
    def items(self):
        """Get the dictionary items.
        
        Returns:
            An iterable of (key, value) pairs
        """
        self._load_dict()
        return self._dict.items()
    
    def __len__(self):
        """Get the number of items in the dictionary.
        
        Returns:
            int: Number of items
        """
        self._load_dict()
        return len(self._dict)


def convert_to_shared_dict(original_dict: Dict[str, Any], dict_id: Optional[str] = None) -> SharedDict:
    """Convert a dictionary to a fast read-only shared memory structure.
    
    This implementation provides much faster read access than traditional 
    shared dictionaries by using direct shared memory access. This approach 
    is optimized for read-heavy workloads where the dictionary doesn't 
    change after creation.
    
    For the best performance (close to local dictionary speed), each worker
    process loads the dictionary into local memory on first access, giving
    near-native speed for all subsequent lookups while still saving memory
    across the main process.
    
    Args:
        original_dict (dict): The regular Python dictionary to convert
        dict_id (str, optional): Unique identifier for this shared dictionary.
                                If None, a UUID will be generated.
        
    Returns:
        SharedDict: A shared memory dictionary with fast read access
    """
    # Generate a unique ID if none provided
    if dict_id is None:
        dict_id = f"shared_dict_{uuid.uuid4().hex}"
    
    print(f"Converting dictionary '{dict_id}' to shared memory...")
    
    # Serialize the dictionary with pickle
    serialized_dict = pickle.dumps(original_dict)
    dict_size_bytes = len(serialized_dict)
    
    # Create a shared memory block to hold the serialized dictionary
    shm = shared_memory.SharedMemory(create=True, size=dict_size_bytes)
    
    # Copy the serialized dictionary to shared memory
    shm.buf[:dict_size_bytes] = serialized_dict
    
    # Store essential information for accessing the shared memory
    shared_dict_info = {
        'name': shm.name,
        'size': dict_size_bytes
    }
    
    # Keep reference to prevent garbage collection
    global _shared_memory_blocks
    _shared_memory_blocks[dict_id] = shm
    
    # Create a wrapper that provides dictionary-like access
    result = SharedDict(shared_dict_info, dict_id)
    
    # Clear the original dictionary to free memory
    original_dict.clear()
    
    # Force garbage collection to release memory
    import gc
    gc.collect()
    
    print(f"Shared memory conversion complete for '{dict_id}' ({dict_size_bytes/1024/1024:.1f} MB)")
    
    return result


def release_shared_dict(dict_id: str) -> bool:
    """Release a specific shared dictionary by its ID.
    
    Args:
        dict_id (str): The ID of the shared dictionary to release
        
    Returns:
        bool: True if the dictionary was successfully released, False otherwise
    """
    global _shared_memory_blocks
    
    if dict_id in _shared_memory_blocks:
        try:
            # Get the shared memory block
            shm = _shared_memory_blocks.pop(dict_id)
            # Close and unlink the shared memory
            shm.close()
            shm.unlink()
            return True
        except Exception as e:
            print(f"Error releasing shared dictionary '{dict_id}': {e}")
            # If an error occurred, put it back in the tracking dictionary
            # so it will be cleaned up on exit
            if dict_id not in _shared_memory_blocks and shm is not None:
                _shared_memory_blocks[dict_id] = shm
    
    return False


def get_shared_dict_ids() -> List[str]:
    """Get a list of all shared dictionary IDs.
    
    Returns:
        List[str]: List of shared dictionary IDs
    """
    global _shared_memory_blocks
    return list(_shared_memory_blocks.keys())
