"use client";
import React, { useMemo, useState } from "react";
import MerckNavbar from "../../_components/MerckNavbar";
import {
  Box,
  Button,
  Card,
  CssBaseline,
  Toolbar,
  Typography,
} from "@mui/material";
import PropertyControls from "@/app/_components/PropertyControls";
import StructureOutput from "@/app/_components/StructureOuptut";
import MolRender from "@/app/_components/MolRender";
import "ketcher-react/dist/index.css";
import dynamic from "next/dynamic";

const propertyNameToKey: { [key: string]: string } = {
  "logP Min": "log_p_min",
  "logP Max": "log_p_max",
  "num molecules": "num_molecules",
  "qed Min": "qed_min",
  "qed Max": "qed_max",
};

const vaeGan: React.FC = () => {
  // eslint-disable-next-line react-hooks/rules-of-hooks
  const [apiResponse, setApiResponse] = useState<any>(null);
  const [genID, setGenID] = useState<any>(null);
  const [inputSmile, setInputSmile] = useState<string>("");

  const KetcherComponent = useMemo(
    () =>
      dynamic(() => import("@/app/_components/KetcherComponent"), {
        ssr: false,
      }),
    []
  );
  const handleGenerateMolecules = async () => {
    const payload = {
      input_smile: inputSmile,
      logP: [parseFloat(propertyValues["logP Min"]), parseFloat(propertyValues["logP Max"])],
      num_molecules: parseFloat(propertyValues["num molecules"]),
      qed: [parseFloat(propertyValues["qed Min"]), parseFloat(propertyValues["qed Max"])],
    };
    //payload["scaffold_smile"] = smile;
    const data = {
      model_name: "vae-gan",
      payload,
    };

    async function polling(pingData: any) {
      const response = await fetch(
        "https://p72f5klivypjgffz23p43th3fy0zweej.lambda-url.us-east-1.on.aws/",
        {
          method: "POST",
          body: JSON.stringify(pingData)
        }
      );
      const responseData = await response.json()

      if (response.ok && responseData.error === "Model output not found. Try polling again after some time.") {
        return polling(pingData);
      } else if (response.ok) {
        setApiResponse(responseData)
        console.log(apiResponse)
      } else {
        console.log("error")
      }
    }


    console.log("Sending payload", data);
    const generationID = await fetch(
      "https://3t777zoaqfoasdu76g335hq37a0uevko.lambda-url.us-east-1.on.aws/",
      {
        method: "POST",
        body: JSON.stringify(data),
      }
    );
    const gen_json = await generationID.json()
    setGenID(gen_json.generation_id)
    console.log(genID);

    const pingData = {
      generation_id: genID
    }

    polling(pingData)
    
  };
    

  // eslint-disable-next-line react-hooks/rules-of-hooks
  const [propertyValues, setPropertyValues] = useState({
    "logP Min": "0",
    "logP Max": "10",
    "qed Min": "0",
    "qed Max": "1",
    "num molecules": "10",
  });

  const aspring = "CC(=O)OC1=CC=CC=C1C(=O)O";

  const handlePropertyChange = (property: string, value: string) => {
    setPropertyValues((prevValues) => ({
      ...prevValues,
      [property]: value,
    }));
  };

  return (
    <Box className="flex">
      <CssBaseline />
      <MerckNavbar />
      <Box component="main" className="flex flex-col flex-grow">
        <Toolbar />
        <Box className="flex flex-row justify-evenly">
          <Box className="w-8/12 p-3" style={{ height: "50vh" }}>
            {" "}
            {/* Keep the existing height */}
            <Box style={{ height: "20px" }}></Box>
            <Button
              variant="outlined"
              onClick={handleGenerateMolecules}
              style={{ marginBottom: "20px" }}
            >
              Generate Molecules
            </Button>
            {apiResponse && (
              <Box className="flex flex-row justify-center flex-wrap">
                <StructureOutput
                  response={apiResponse}
                  isMultiObj={false}
                  input_smile={inputSmile}
                />
              </Box>
            )}
          </Box>
          <Card className="w-3/12 p-3 overflow-scroll">
            <PropertyControls
              propertyValues={propertyValues}
              handlePropertyChange={handlePropertyChange}
            />
          </Card>
        </Box>
      </Box>
    </Box>
  );
};

export default vaeGan;
