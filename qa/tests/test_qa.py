"""
Unit and regression test for the qa package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import qa


def test_qa_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "qa" in sys.modules
