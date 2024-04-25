from uuid import uuid4
import boto3
from botocore.exceptions import ClientError
import json
import logging

MODEL_LIST = ["vae-gan", "multiobj-rationale", "scaffold-constrained"]
TABLE_NAME = "MoleculeGenerationStorage"
QUEUE_NAME = "generation_queue.fifo"

sqs = boto3.resource("sqs")
dynamodb = boto3.client("dynamodb")
queue = sqs.get_queue_by_name(QueueName=QUEUE_NAME)
logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def lambda_handler(event, context):
    logger.info("Received event: " + json.dumps(event, indent=2))

    data = json.loads(json.dumps(event))
    data = json.loads(data["body"])
    if "payload" not in data:
        return {
            "statusCode": 400,
            "body": json.dumps("error: Missing payload in request."),
        }
    payload = data["payload"]
    if "model_name" not in data:
        return {
            "statusCode": 400,
            "body": json.dumps("error: Missing model name in request."),
        }
    model_name = data["model_name"]
    if model_name not in MODEL_LIST:
        return {
            "statusCode": 400,
            "body": json.dumps("error: Invalid model name specified."),
        }
    try:
        generation_id = uuid4().hex
        response = queue.send_message(
            MessageBody=json.dumps(data),
            MessageAttributes={
                "generation_id": {"StringValue": generation_id, "DataType": "String"}
            },
            MessageGroupId=model_name,
        )
    except ClientError as error:
        logger.exception("Send message failed: %s", payload)
        raise error
    else:
        dynamodb.put_item(
            TableName=TABLE_NAME,
            Item={
                "generation_id": {"S": generation_id},
            },
        )
        return {
            "statusCode": 200,
            "body": json.dumps(
                {
                    "generation_id": generation_id,
                    "status": "Generation queued successfully.",
                }
            ),
        }
