import subprocess
from subprocess import CalledProcessError
import sagemaker_model_handler
from sagemaker_inference import model_server
import os


def _start_mms():
    # by default the number of workers per model is 1, but we can configure it through the
    # environment variable below if desired.
    # os.environ['SAGEMAKER_MODEL_SERVER_WORKERS'] = '2'
    print("Starting MMS -> running ", sagemaker_model_handler.__file__)
    model_server.start_model_server(
        handler_service=sagemaker_model_handler.__file__ + ":handle"
    )


def main():
    _start_mms()
    # prevent docker exit
    subprocess.call(["tail", "-f", "/dev/null"])


main()
