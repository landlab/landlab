def pytest_generate_tests(metafunc):
    if "format" in metafunc.fixturenames:
        metafunc.parametrize(
            "format", ("NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_CLASSIC", "NETCDF3_64BIT")
        )
    elif "at" in metafunc.fixturenames:
        metafunc.parametrize("at", ("node", "link", "patch", "corner", "face", "cell"))
