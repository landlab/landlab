#! /usr/bin/env python
import os

import pytest

from landlab.core.messages import assert_or_print
from landlab.core.messages import error_message
from landlab.core.messages import format_message
from landlab.core.messages import split_paragraphs
from landlab.core.messages import warning_message

LOREM_IPSUM = os.linesep.join(
    [
        "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.",  # noqa: B950
        "",
        "Pharetra pharetra massa massa ultricies mi quis hendrerit.",
        "",
        "Dictumst vestibulum rhoncus est pellentesque. Sed viverra tellus in hac habitasse platea dictumst vestibulum rhoncus.",  # noqa: B950
    ]
)


def test_split_paragraphs_cr():
    """Test splitting paragraphs with carriage returns."""
    text = """
Pharetra pharetra massa massa ultricies mi quis hendrerit.\r\rDictumst vestibulum rhoncus est pellentesque.
    """  # noqa: B950
    assert split_paragraphs(text, linesep="\r") == [
        "Pharetra pharetra massa massa ultricies mi quis hendrerit.",
        "Dictumst vestibulum rhoncus est pellentesque.",
    ]


def test_split_paragraphs_lf():
    """Test splitting paragraphs with line feeds."""
    text = """
Pharetra pharetra massa massa ultricies mi quis hendrerit.\n\nDictumst vestibulum rhoncus est pellentesque.
    """  # noqa: B950
    assert split_paragraphs(text, linesep="\n") == [
        "Pharetra pharetra massa massa ultricies mi quis hendrerit.",
        "Dictumst vestibulum rhoncus est pellentesque.",
    ]


def test_split_paragraphs_crlf():
    """Test splitting paragraphs with carriage returns and line feeds."""
    text = """
Pharetra pharetra massa massa ultricies mi quis hendrerit.\r\n\r\nDictumst vestibulum rhoncus est pellentesque.
    """  # noqa: B950
    assert split_paragraphs(text, linesep="\r\n") == [
        "Pharetra pharetra massa massa ultricies mi quis hendrerit.",
        "Dictumst vestibulum rhoncus est pellentesque.",
    ]


def test_empty_message():
    """Test formatting an empty string."""
    assert format_message("") == ""


def test_one_line():
    """Test a single line message."""
    assert format_message("lorem ipsum") == "lorem ipsum"


def test_leading_whitespace():
    """Test a single line message."""
    assert format_message("   lorem ipsum") == "lorem ipsum"


def test_one_long_line():
    """Test a line that needs to be wrapped."""
    msg = (
        "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod "
        "tempor incididunt ut labore et dolore magna aliqua."
    )

    assert format_message(msg) == os.linesep.join(
        [
            "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do",
            "eiusmod tempor incididunt ut labore et dolore magna aliqua.",
        ]
    )


def test_multiline():
    msg = """lorem
ipsum
    """
    assert format_message(msg) == "lorem ipsum"


def test_multiple_paragraphs():
    assert format_message(LOREM_IPSUM) == os.linesep.join(
        [
            "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do",
            "eiusmod tempor incididunt ut labore et dolore magna aliqua.",
            "",
            "Pharetra pharetra massa massa ultricies mi quis hendrerit.",
            "",
            "Dictumst vestibulum rhoncus est pellentesque. Sed viverra tellus in",
            "hac habitasse platea dictumst vestibulum rhoncus.",
        ]
    )


def test_warning_message():
    msg = "Pharetra pharetra massa massa ultricies mi quis hendrerit."

    assert warning_message(msg) == os.linesep.join(
        [
            "WARNING",
            "=======",
            "",
            "Pharetra pharetra massa massa ultricies mi quis hendrerit.",
        ]
    )


def test_error_message():
    msg = "Pharetra pharetra massa massa ultricies mi quis hendrerit."

    assert error_message(msg) == os.linesep.join(
        [
            "ERROR",
            "=====",
            "",
            "Pharetra pharetra massa massa ultricies mi quis hendrerit.",
        ]
    )


def test_warning_message_is_none():
    assert warning_message() == os.linesep.join(["WARNING", "======="])


def test_error_message_is_none():
    assert error_message() == os.linesep.join(["ERROR", "====="])


def test_assert_or_pass():
    assert_or_print(True, onerror="pass")
    assert_or_print(False, onerror="pass")


def test_assert_or_warn():
    assert_or_print(True, onerror="warn")
    assert_or_print(False, onerror="warn")


def test_assert_or_raise():
    assert_or_print(True, onerror="raise")
    with pytest.raises(AssertionError):
        assert_or_print(False, onerror="raise")
