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
      log_p_min: parseFloat(propertyValues["logP Min"]),
      log_p_max: parseFloat(propertyValues["logP Max"]),
      num_molecules: parseFloat(propertyValues["num molecules"]),
      qed_min: parseFloat(propertyValues["qed Min"]),
      qed_max: parseFloat(propertyValues["qed Max"]),
    };
    //payload["scaffold_smile"] = smile;
    const data = {
      model_type: "vae-gan",
      payload,
    };

    console.log("Sending payload", data);
    const response = await fetch(
      "https://ezu74lbfo2imcxnmbblg3hkhqq0oqtxo.lambda-url.us-east-1.on.aws/",
      {
        method: "POST",
        body: JSON.stringify(data),
      }
    );
    setApiResponse(await response.json());
    console.log(apiResponse);
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
                  response={apiResponse.filtered_smiles}
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
