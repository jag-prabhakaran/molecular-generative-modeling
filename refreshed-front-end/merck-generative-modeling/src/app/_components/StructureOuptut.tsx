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
  
  

  const StructureOutput: React.FC<Molecule[]> = ({ respone }) => {
    {respone.map((molecule: object) => {
        return (
            <MolRender
            key={molecule.smile}
            molecule={molecule} />
        )
    })}
  };


export default StructureOutput