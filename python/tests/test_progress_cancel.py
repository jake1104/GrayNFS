import graynfs


def test_progress_false_cancels() -> None:
    seen: list[dict] = []

    def cancel(event: dict) -> bool:
        seen.append(event)
        return False

    try:
        graynfs.factor_with_progress(15, cancel)
    except graynfs.GrayNFSError as exc:
        assert exc.category == "cancelled"
        assert exc.detail == "the operation was cancelled"
        assert str(exc) == "cancelled: the operation was cancelled"
    else:
        raise AssertionError("expected GrayNFSError")
    assert seen


def test_progress_true_continues() -> None:
    assert graynfs.factor_with_progress(15, lambda event: True) == [3, 5]


def test_progress_none_continues() -> None:
    assert graynfs.factor_with_progress(15, lambda event: None) == [3, 5]
