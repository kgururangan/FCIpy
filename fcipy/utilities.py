import os
import psutil

def get_memory_usage():
    """Displays the percentage of used RAM and available memory. Useful for
    investigating the memory usages of various routines.
    Usage:
    >> # gives a single float value
    >> psutil.cpu_percent()
    >> # gives an object with many fields
    >> psutil.virtual_memory()
    >> # you can convert that object to a dictionary
    >> dict(psutil.virtual_memory()._asdict())
    >> # you can have the percentage of used RAM
    >> psutil.virtual_memory().percent
    >> 79.2
    >> # you can calculate percentage of available memory
    >> psutil.virtual_memory().available * 100 / psutil.virtual_memory().total
    >> 20.
    """
    current_process = psutil.Process(os.getpid())
    memory = current_process.memory_info().rss # RSS (e.g., RAM usage) memory in bytes
    return memory / (1024 * 1024)
