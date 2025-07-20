import importlib

def test_import_cvector():
    from cmtj import CVector
    assert CVector.__name__ == "CVector"

def test_import_fieldscan():
    from cmtj.utils import FieldScan
    assert hasattr(FieldScan, "angle2vector")
