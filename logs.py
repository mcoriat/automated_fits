import logging
import sys
import traceback
import warnings

from pathlib import Path


class UltranestFilter(logging.Filter):
    def filter(self, record):
        return not record.getMessage().startswith("iteration=")


class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, level):
       self.logger = logger
       self.level = level
       self.linebuf = ''

    def write(self, buf):
       for line in buf.rstrip().splitlines():
          self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass


def set_logger(srcid, model, stdout_to_log=True, fmt=None):
    log_file = Path("logs", model, f"fit_{srcid}.log")
    filehandler = logging.FileHandler(log_file, "w")                                                                                                                                                       
                                                                                                                                                                                                           
    if fmt is None:                                                                                                                                                                                        
        #infofmt = "%(levelname)s:%(asctime)s: %(module)s:%(funcName)s: %(message)s"
        infofmt = "[%(name)s %(levelname)s]: %(message)s"
        fmt = logging.Formatter(infofmt, datefmt="%I:%M:%S")

    filehandler.setFormatter(fmt)
    filehandler.addFilter(UltranestFilter())

    # root logger - Good to get it only once.
    logger = logging.getLogger()

    # remove the existing file handlers
    for hdlr in logger.handlers[:]:
        logger.removeHandler(hdlr)

    logger.addHandler(filehandler)
    logger.setLevel(logging.INFO)
    
    if stdout_to_log:
        sys.stdout = StreamToLogger(logger, logging.INFO)
        sys.stderr = StreamToLogger(logger, logging.ERROR)

    logger = logging.getLogger("sherpa")
    logger.setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", message='displayed errorbars')
    
    logger = logging.getLogger("ultranest")
    logger.setLevel(logging.INFO)


def log_exception(exception):
    logging.error(''.join(traceback.format_tb(exception.__traceback__)))
    logging.error(exception)
