import enum

class Precision(enum.Enum):
    DEFAULT = 0
    TF32 = 1
    TF32x3 = 2
    IEEE = 3
