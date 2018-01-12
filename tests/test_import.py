import fair
version = fair.__version__

# Are we dealing with v1.1-x?
def test_import():
    assert version[:3] == '1.1'
