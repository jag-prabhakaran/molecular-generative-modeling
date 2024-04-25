import boto3
from botocore.exceptions import ClientError
import json
import logging

TABLE_NAME = "MoleculeGenerationStorage"

dynamodb = boto3.client("dynamodb")
logging.basicConfig()
logger = logging.getLogger()
logger.setLevel(logging.INFO)


def lambda_handler(event, context):
    logger.info("Received event: " + json.dumps(event, indent=2))

    data = json.loads(json.dumps(event))
    data = json.loads(data["body"])

    if "generation_id" not in data:
        return {
            "statusCode": 400,
            "body": json.dumps("error: Missing generation_id in request."),
        }
    generation_id = data["generation_id"]
    try:
        molecule_row = dynamodb.get_item(
            TableName=TABLE_NAME,
            Key={"generation_id": {"S": generation_id}},
        )
        if "model_output" not in molecule_row["Item"]:
            return {
                "statusCode": 500,
                "body": json.dumps(
                    {
                        "error": "Model output not found. Try polling again after some time."
                    }
                ),
            }
        model_json = molecule_row["Item"]["model_output"]["S"]
    except ClientError as error:
        logger.exception("Get item failed: %s", generation_id)
        raise error
    else:
        return {
            "statusCode": 200,
            "body": model_json,
        }
