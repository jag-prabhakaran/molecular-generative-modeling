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
  const { key, molecule } = props;
  return (
    <Box>
      <MoleculeStructure structure={molecule.smile} id={molecule.smile} />
    </Box>
  );
};

export default MolRender;
