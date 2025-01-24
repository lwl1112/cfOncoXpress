# Marjorie Roskes (mcl4001@med.cornell.edu)

import logging
import os

logformat = "%(message)s at time %(asctime)s."

subjob = os.environ.get("LOG_SUBJOB", None)

if subjob is not None:
    logformat = f"Subjob {subjob}: {logformat}"

logging.basicConfig(format=logformat, datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger("Logger")

logger.setLevel(logging.DEBUG)