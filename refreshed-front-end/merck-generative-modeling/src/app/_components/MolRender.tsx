"use client";
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
import React, { HTMLAttributes, useState } from "react";
import MoleculeStructure from "./MoleculeStructure";
import Molecule from "./MoleculeType";

const MolRender = (props: { key: string; molecule: Molecule }) => {
  const round = (value: number) => Math.round(value * 1000) / 1000;
  const { key, molecule } = props;
  return (
    <Box>
      <MoleculeStructure
        structure={molecule.smile}
        id={molecule.smile}
        svgMode={true}
      />
      <Typography>LogP: {round(molecule.logP)}</Typography>
      <Typography>QED score: {round(molecule.qed)}</Typography>
      <Typography>Molecular Weight: {round(molecule.mol_weight)}</Typography>
      <Typography>
        Number of H donors: {round(molecule.num_h_donors)}
      </Typography>
    </Box>
  );
};

export default MolRender;
