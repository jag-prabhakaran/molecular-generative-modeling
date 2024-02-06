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

  interface Molecule {
    "smile": string,
    "logP": number,
    "qed": number,
    "mol_weight": number,
    "num_h_donors": number
  }
  
import MolRender from './MolRender';

const StructureOutput: React.FC<{ response: Molecule[] }> = ({response}) => {
  return (
    <>
      {response.map((molecule: Molecule) => {
        return (
          <MolRender
            key={molecule.smile}
            molecule={molecule}
          />
        );
      })}
    </>
  );
};


export default StructureOutput