"use client";
import React, { use, useMemo, useState } from "react";
import MerckNavbar from "../../_components/MerckNavbar";
import {
  Box,
  Button,
  Card,
  CssBaseline,
  Toolbar,
  Typography,
} from "@mui/material";
import dynamic from "next/dynamic";
import PropertyControls from "@/app/_components/PropertyControls";
import StructureOutput from "@/app/_components/StructureOuptut";
import MolRender from "@/app/_components/MolRender";
import Molecule from "@/app/_components/MoleculeType";

const propertyNameToKey: { [key: string]: string } = {
  "logP Min": "log_p_min",
  "logP Max": "log_p_max",
  "upper bound": "num_molecules",
  "qed Min": "qed_min",
  "qed Max": "qed_max",
};



const vaeGan: React.FC = () => {
  // eslint-disable-next-line react-hooks/rules-of-hooks
  const KetcherComponent = useMemo(
    () =>
      dynamic(() => import("@/app/_components/KetcherComponent"), {
        ssr: false,
      }),
    []
  );
  const [apiResponse, setApiResponse] = useState<any>(null);
  const [inputSmile, setInputSmile] = useState<string>("");

  const handleGenerateMolecules = async () => {
    const smile = await (window as any).ketcher.getSmiles();
    setInputSmile(smile.replaceAll(":", ""));

    console.log(inputSmile);
    const payload = {
      log_p_min: parseFloat(propertyValues["logP Min"]),
      log_p_max: parseFloat(propertyValues["logP Max"]),
      num_molecules: parseFloat(propertyValues["upper bound"]),
      qed_min: parseFloat(propertyValues["qed Min"]),
      qed_max: parseFloat(propertyValues["qed Max"]),
      scaffold_smile: smile.replaceAll(":", ""),
    };

    function removeDuplicateMols(array: [Molecule]) {
      const uniqueObjects = new Set(
        array.map((obj: any) => JSON.stringify(obj))
      );
      return Array.from(uniqueObjects).map((str) => JSON.parse(str));
    }

    const data = {
      model_type: "scaffold-constrained",
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
    const res_json = await response.json();
    res_json.filtered_smiles = removeDuplicateMols(res_json.filtered_smiles);
    setApiResponse(res_json);
    console.log(apiResponse);
  };

  // eslint-disable-next-line react-hooks/rules-of-hooks
  const [propertyValues, setPropertyValues] = useState({
    "logP Min": "0",
    "logP Max": "10",
    "qed Min": "0",
    "qed Max": "1",
    "upper bound": "10",
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
            <KetcherComponent />
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
