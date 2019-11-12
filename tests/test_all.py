import testing_classes as tc
from testing_classes import testing_cleaning

testing_cleaning()  # This function is not a real test. Is just to clean the testing area before and after testing


def test_1():
    tc.Convert2SDFTest().test_files_creation()


def test_2():
    tc.Convert2SDFTest().test_file_content()


def test_3():
    tc.Prepare2FragTest().test_file_creation()


def test_4():
    tc.Prepare2FragTest().test_file_content()


def test_5():
    tc.DetectorTest().test_tree()


testing_cleaning()

