"""Thread-safe rate limiter and exponential backoff retry mechanisms."""

import functools
import random
import time
from typing import Any, Callable, TypeVar

from shared.logging.logger import get_logger

logger = get_logger("shared.rate_limiter")
F = TypeVar("F", bound=Callable[..., Any])


class RateLimiter:
    """Token-bucket or minimal interval rate limiter for polite API access."""

    def __init__(self, max_calls_per_second: float = 3.0):
        self.min_interval = 1.0 / max_calls_per_second
        self.last_call_time = 0.0

    def wait(self) -> None:
        """Sleep if necessary to respect the minimal interval between calls."""
        now = time.time()
        elapsed = now - self.last_call_time
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self.last_call_time = time.time()


def retry_with_backoff(
    max_retries: int = 4,
    base_delay: float = 2.0,
    max_delay: float = 30.0,
    jitter: bool = True,
    exceptions: tuple = (Exception,),
) -> Callable[[F], F]:
    """Decorator for retrying unstable network functions with exponential backoff and jitter."""

    def decorator(func: F) -> F:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            attempt = 0
            while True:
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    attempt += 1
                    if attempt > max_retries:
                        logger.error(
                            f"Function {func.__name__} failed after {max_retries} retries: {e}"
                        )
                        raise
                    delay = min(max_delay, base_delay * (2 ** (attempt - 1)))
                    if jitter:
                        delay += random.uniform(0.1, 0.5 * delay)
                    logger.warning(
                        f"Attempt {attempt}/{max_retries} for {func.__name__} failed ({e}). Retrying in {delay:.2f}s..."
                    )
                    time.sleep(delay)

        return wrapper  # type: ignore

    return decorator
