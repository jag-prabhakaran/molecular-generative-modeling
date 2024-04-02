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
import initRDKitModule from "../../../../public/js/RDKit_minimal";
import { RDKitLoader } from "@rdkit/rdkit";
import "ketcher-react/dist/index.css";
import dynamic from "next/dynamic";

const propertyNameToKey: { [key: string]: string } = {
  "logP Min": "log_p_min",
  "logP Max": "log_p_max",
  "upper bound": "num_molecules",
  "qed Min": "qed_min",
  "qed Max": "qed_max",
};

const vaeGan: React.FC = () => {
  const KetcherComponent = useMemo(
    () =>
      dynamic(() => import("@/app/_components/KetcherComponent"), {
        ssr: false,
      }),
    []
  );
  // eslint-disable-next-line react-hooks/rules-of-hooks
  const [apiResponse, setApiResponse] = useState<any>(null);

  const initRDKit = (() => {
    let rdkitLoadingPromise: Promise<any> | undefined = undefined;

    return () => {
      /**
       * Utility function ensuring there's only one call made to load RDKit
       * It returns a promise with the resolved RDKit API as value on success,
       * and a rejected promise with the error on failure.
       *
       * The RDKit API is also attached to the global object on successful load.
       */
      if (!rdkitLoadingPromise) {
        rdkitLoadingPromise = new Promise((resolve, reject) => {
          initRDKitModule({ locateFile: (file: string) => `/js/${file}` })
            .then((RDKit: RDKitLoader) => {
              resolve(RDKit);
            })
            .catch((e: Error) => {
              reject();
            });
        });
      }

      return rdkitLoadingPromise;
    };
  })();
  var RDKitModule: any;
  initRDKit().then((RDKit: any) => {
    RDKitModule = RDKit;
  });

  const [inputSmile, setInputSmile] = useState<string>("");

  const handleGenerateMolecules = async () => {
    let smile: string = await (window as any).ketcher.getSmiles();
    smile = RDKitModule.get_mol(smile).get_smiles();
    const regex = /([a-zA-z][0-9]*)(\*)/gm;
    smile = smile.replaceAll(regex, "[$1:1]");
    console.log(smile);
    setInputSmile(smile)

    const payload = {
      log_p_min: parseFloat(propertyValues["logP Min"]),
      log_p_max: parseFloat(propertyValues["logP Max"]),
      num_molecules: parseFloat(propertyValues["num molecules"]),
      qed_min: parseFloat(propertyValues["qed Min"]),
      qed_max: parseFloat(propertyValues["qed Max"]),
      rationale: [smile],
    };

    //payload["scaffold_smile"] = smile;
    const data = {
      model_type: "multiobj-rationale",
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
    setApiResponse(res_json);
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
                  response={apiResponse.filtered_output_objects}
                  isMultiObj={true}
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
