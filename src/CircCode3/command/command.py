# -*- coding = utf-8 -*-
import subprocess
from ..utils.logs import logger


class Command:
    def __init__(self, software: str, *args):
        self.command = [str(i) for i in args]
        self.software = software

    def to_command(self):
        command = self.software + " " + " ".join(self.command)
        return command

    def command_running(self):
        """ running software build index
        """
        logger.info(f"Start running command:\n{self.to_command()}")
        result = subprocess.run(self.to_command(), capture_output=True, check=True, shell=True)
        if not result.returncode:
            logger.info(f"{self.software} running completed!")
        else:
            logger.error(f"{self.software} running error:\n{result.stderr.decode()}")
            logger.error(f"{result.stdout.decode()}")
            exit(1)
