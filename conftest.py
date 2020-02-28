def pytest_addoption(parser):
    parser.addoption(
        "--run-notebook", action="store_true", default=False, help="run notebook tests"
    )


