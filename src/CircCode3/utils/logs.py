# -*- coding = utf-8 -*-
import logging
import sys


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers= [
        logging.FileHandler("./CircCode3.log", mode="w"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger("CircCode3")

def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error("ERROR:", exc_info=(exc_type, exc_value, exc_traceback))
    sys.exit(1)

sys.excepthook = handle_exception

