
def test_import_cvector():
    from cmtj import CVector
    assert CVector.__name__ == "CVector"
    from cmtj.stack import ParallelStack, SeriesStack
    assert hasattr(ParallelStack, "runSimulation")
    assert hasattr(SeriesStack, "runSimulation")
    from cmtj.reservoir import GroupInteraction
    assert hasattr(GroupInteraction, "runSimulation")
    from cmtj.llgb import LLGBJunction
    assert hasattr(LLGBJunction, "runSimulation")

def test_import_fieldscan():
    from cmtj.utils import FieldScan
    assert hasattr(FieldScan, "angle2vector")