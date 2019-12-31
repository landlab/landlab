def pytest_generate_tests(metafunc):
    if "at" in metafunc.fixturenames:
        metafunc.parametrize("at", ("node", "link", "patch", "corner", "face", "cell"))
