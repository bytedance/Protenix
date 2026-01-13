import requests
import json
import os
import time

# Configuration
API_KEY = os.environ.get("RUNPOD_API_KEY", "YOUR_API_KEY")
ENDPOINT_ID = os.environ.get("RUNPOD_ENDPOINT_ID", "YOUR_ENDPOINT_ID")

# use /runsync for synchronous or /run for async
URL = f"https://api.runpod.ai/v2/{ENDPOINT_ID}/runsync"

headers = {
    "Content-Type": "application/json",
    "Authorization": f"Bearer {API_KEY}"
}

# Payload matching Protenix input format
data = {
    "input": {
        "name": "test_request_py",
        "sequences": [
            {
                "proteinChain": {
                    "sequence": "MAEVIRSSAFWRSFPIFEEFDSETLCELSGIASYRKWSAGTVIFQRGDQGDYMIVVVSGRIKLSLFTPQGRELMLRQHEAGALFGEMALLDGQPRSADATAVTAAEGYVIGKKDFLALITQRPKTAEAVIRFLCAQLRDTTDRLETIALYDLNARVARFFLATLRQIHGSEMPQSANLRLTLSQTDIASILGASRPKVNRAILSLEESGAIKRADGIICCNVGRLLSIADPEEDLEHHHHHHHH",
                    "count": 1
                }
            }
        ],
        "use_msa": False
    }
}

print(f"Sending request to {URL}...")
try:
    response = requests.post(URL, headers=headers, json=data)
    response.raise_for_status()
    
    result = response.json()
    print("Response Status Code:", response.status_code)
    
    if "output" in result:
        print("Job Output Received!")
        # If output_files are present, maybe save them?
        output_data = result.get("output", {})
        if "output_files" in output_data:
             print(f"Received {len(output_data['output_files'])} file(s).")
             for fname, content in output_data["output_files"].items():
                 print(f"Saving {fname}...")
                 with open(fname, "w") as f:
                     f.write(content)
        else:
            print("Output:", json.dumps(output_data, indent=2))
    else:
        print("Full Response:", json.dumps(result, indent=2))

except requests.exceptions.RequestException as e:
    print(f"Error: {e}")
    if hasattr(e, 'response') and e.response is not None:
        print("Response text:", e.response.text)
