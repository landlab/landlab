from collections.abc import Mapping
from typing import Any

from landlab.core.messages import deprecation_message


class DeprecatedSignature(DeprecationWarning):
    def __init__(
        self,
        name: str,
        new: tuple[tuple[Any, ...], Mapping[str, Any]] | None = None,
    ):
        super().__init__(name, new)

    def __str__(self) -> str:
        name, new = self.args

        if new is None:
            use = None
        else:
            new = self._construct_call(name, new[0], new[1])
            use = f">>> grid = {new}"

        return deprecation_message(
            "You are using a deprecated calling signature.", use=use
        )

    @staticmethod
    def _construct_call(
        name: str,
        args: tuple[Any, ...],
        kwds: Mapping[str, Any],
    ) -> str:
        signature = ", ".join(
            [repr(arg) for arg in args] + [f"{k}={repr(v)}" for k, v in kwds.items()]
        )
        return f"{name}({signature})"
