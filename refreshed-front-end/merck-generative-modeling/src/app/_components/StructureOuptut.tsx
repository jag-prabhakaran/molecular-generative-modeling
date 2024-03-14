import {
  Box,
  Button,
  Card,
  CardActions,
  CardContent,
  Input,
  List,
  Slider,
  TextField,
  Typography,
} from "@mui/material";
import React, { useState } from "react";
import Molecule from "./MoleculeType";

import MolRender from "./MolRender";

const StructureOutput: React.FC<{
  response: Molecule[];
  isMultiObj: boolean;
}> = ({ response, isMultiObj }) => {
  return (
    <>
      {isMultiObj
        ? response.map((molecule: Molecule) => {
            return (
              <MolRender
                key={molecule.output_smile}
                molecule={molecule}
                isMultiObj={true}
                smile={molecule.output_smile}
              />
            );
          })
        : response.map((molecule: Molecule) => {
            return (
              <MolRender
                key={molecule.smile}
                molecule={molecule}
                isMultiObj={false}
                smile={molecule.smile}
              />
            );
          })}
    </>
  );
};

export default StructureOutput;
