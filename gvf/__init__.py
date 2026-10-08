"""gvf - Gradually Varied Flow and hydraulic jump profiler for open channels."""
from .profile import compute
from .sections import ChannelError

__all__ = ["compute", "ChannelError"]
__version__ = "1.0.0"
