import os
import pkgutil


def test_import_main():
    """
    This imports PyThea
    """
    no_requests = False
    try:
        pass
    except ImportError:
        no_requests = True
    assert no_requests is False


def test_imports_all():
    """
    This imports all modules in PyThea
    """
    def on_error():
        try:
            raise
        except Warning:
            pass

    test_directory = os.path.dirname(__file__)

    for imper, nm, ispkg in pkgutil.walk_packages([f'{test_directory}'], 'PyThea.',
                                                  onerror=on_error):
        imper.find_spec(nm)
