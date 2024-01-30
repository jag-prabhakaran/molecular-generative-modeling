"use client";
import React, { useState } from "react";
import MerckNavbar from "../../_components/MerckNavbar";
import {
  Box,
  Button,
  Card,
  CssBaseline,
  Toolbar,
  Typography,
} from "@mui/material";
import KetcherComponent from "@/app/_components/KetcherComponent";
import PropertyControls from "@/app/_components/PropertyControls";

const propertyNameToKey: { [key: string]: string } = {
  "logP Min": "log_p_min",
  "logP Max": "log_p_max",
  "num molecules": "num_molecules",
//   "qed Min": "qed_min",
//   "qed Max": "qed_max",
};

const vaeGan: React.FC = () => {


    // eslint-disable-next-line react-hooks/rules-of-hooks
    const [apiResponse, setApiResponse] = useState<any>(null);

    const handleGenerateMolecules = async () => {
      const payload = Object.keys(propertyValues).reduce(
        (result: { [key: string]: any }, property) => {
          const key = propertyNameToKey[property];
          if (key) {
            result[key] = propertyValues[property];
          }
          return result;
        },
        {}
      );

      const smile = await (window as any).ketcher.getSmiles();
    //   payload["scaffold_smile"] = smile;
      const data = {
        model_type: "vae-gan",
        payload,
      };

      console.log("Sending payload", data);
      const response = await (
        await fetch(
          "https://ezu74lbfo2imcxnmbblg3hkhqq0oqtxo.lambda-url.us-east-1.on.aws/",
          {
            method: "POST",
            body: JSON.stringify(data),
          }
        )
      );
      setApiResponse(response.text());
        console.log(apiResponse);
    };

  // eslint-disable-next-line react-hooks/rules-of-hooks
  const [propertyValues, setPropertyValues] = useState<{
    [key: string]: string;
  }>({});

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
            <KetcherComponent />
            {/* <Box className='w-full h-full' bgcolor="red"></Box> */}
          </Box>
          <Card className="w-3/12 p-3 overflow-scroll">
            <PropertyControls
              propertyValues={propertyValues}
              handlePropertyChange={handlePropertyChange}
            />
          </Card>
        </Box>
        <Box className="flex flex-row justify-center">
          <Button variant="contained" onClick={handleGenerateMolecules}>
            Generate Molecules
          </Button>
        </Box>
      </Box>
    </Box>
  );
};

export default vaeGan;
