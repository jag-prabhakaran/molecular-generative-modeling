import boto3
import json

runtime = boto3.client("runtime.sagemaker")  # type: ignore


def lambda_handler(event, context):
    print("Received event: " + json.dumps(event, indent=2))

    data = json.loads(json.dumps(event))
    data = json.loads(data["body"])
    print(data)
    if "payload" not in data:
        return json.dumps({"error": "Missing payload in request."})
    payload = data["payload"]
    print(payload)
    model_type = data["model_type"]
    if model_type == "vae-gan":
        if (
            "num_molecules" not in payload
            or "log_p_min" not in payload
            or "log_p_max" not in payload
            or "num_molecules" not in payload
        ):
            return json.dumps({"error": "Missing parameters in payload."})
    elif model_type == "scaffold-constrained":
        if (
            "scaffold_smile" not in payload
            or "log_p_min" not in payload
            or "log_p_max" not in payload
        ):
            return json.dumps({"error": "Missing parameters in payload."})
    response = runtime.invoke_endpoint(
        EndpointName=model_type,
        ContentType="application/json",
        Body=json.dumps(payload),
    )
    print(response)
    out_json = response["Body"].read().decode("utf-8")
    print(out_json)
    return out_json
