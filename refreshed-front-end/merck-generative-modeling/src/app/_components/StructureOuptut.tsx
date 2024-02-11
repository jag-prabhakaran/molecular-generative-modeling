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

const StructureOutput: React.FC<{ response: Molecule[] }> = ({ response }) => {
  return (
    <>
      {response.map((molecule: Molecule) => {
        return <MolRender key={molecule.smile} molecule={molecule} />;
      })}
    </>
  );
};

export default StructureOutput;
