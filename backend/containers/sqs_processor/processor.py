import logging
import json
from math import log
import boto3
from mypy_boto3_sqs.service_resource import Message
import sys

QUEUE_URL = "https://sqs.us-east-1.amazonaws.com/713411356940/generation_queue.fifo"
TABLE_NAME = "MoleculeGenerationStorage"

sagemaker_runtime = boto3.client("sagemaker-runtime")
sqs = boto3.resource("sqs")
dynamodb = boto3.client("dynamodb")
client = boto3.client("sqs")
queue = sqs.get_queue_by_name(QueueName="generation_queue.fifo")
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler(sys.stdout))


def process_message(message: Message):
    logger.info(f"Got Message: {message.body}")
    message_data = json.loads(message.body)
    generation_id = message.message_attributes.get("generation_id")
    if generation_id:
        model_response = sagemaker_runtime.invoke_endpoint(
            EndpointName=message_data.get("model_name"),
            ContentType="application/json",
            Body=json.dumps(message_data["payload"]),
        )
        model_output_json = model_response["Body"].read().decode("utf-8")
        dynamodb.put_item(
            TableName=TABLE_NAME,
            Item={
                "generation_id": {"S": generation_id.get("StringValue")},
                "model_output": {"S": model_output_json},
            },
        )
        client.delete_message(QueueUrl=QUEUE_URL, ReceiptHandle=message.receipt_handle)


if __name__ == "__main__":

    while True:
        messages = queue.receive_messages(
            MaxNumberOfMessages=1, WaitTimeSeconds=5, MessageAttributeNames=["All"]
        )
        if len(messages) > 0:
            process_message(messages[0])
