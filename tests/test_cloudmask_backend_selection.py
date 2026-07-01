import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import run_fylat


def main() -> None:
    original_env = os.environ.get("FYLAT_CLOUDMASK_BACKEND")
    try:
        os.environ.pop("FYLAT_CLOUDMASK_BACKEND", None)
        assert run_fylat.cloudmask_backend_name() == "fortran"

        os.environ["FYLAT_CLOUDMASK_BACKEND"] = "cpp_ocean_day"
        assert run_fylat.cloudmask_backend_name() == "cpp_ocean_day"

        os.environ["FYLAT_CLOUDMASK_BACKEND"] = "auto"
        assert run_fylat.cloudmask_backend_name() == "auto"

        os.environ["FYLAT_CLOUDMASK_BACKEND"] = "bad"
        try:
            run_fylat.cloudmask_backend_name()
        except ValueError as exc:
            assert "Unsupported FYLAT_CLOUDMASK_BACKEND" in str(exc)
        else:
            raise AssertionError("unsupported cloudmask backend did not raise")
    finally:
        if original_env is None:
            os.environ.pop("FYLAT_CLOUDMASK_BACKEND", None)
        else:
            os.environ["FYLAT_CLOUDMASK_BACKEND"] = original_env


if __name__ == "__main__":
    main()
