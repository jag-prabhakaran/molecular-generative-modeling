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
import { DataGrid, GridColDef, GridRenderCellParams, GridToolbar } from "@mui/x-data-grid";
import MolRender from "./MolRender";
import MoleculeStructure from "./MoleculeStructure";

const StructureOutput: React.FC<{
  response: Molecule[];
  isMultiObj: boolean;
  input_smile: string
}> = ({ response, isMultiObj, input_smile }) => {
  const round = (value: number) => Math.round(value * 1000) / 1000;

  const rows = response.map((molecule, index) => ({
    id: index,
    logP: round(molecule.logP),
    qed: round(molecule.qed),
    mol_weight: round(molecule.mol_weight),
    num_h_donors: round(molecule.num_h_donors),
    moleculeStructure: (
      <MoleculeStructure
        structure={molecule.smile}
        id={`mol_${index}`}
        subStructure={input_smile}
        svgMode={true}
      />
    ),
  }));

  const columns = [
    { field: 'moleculeStructure', headerName: 'Structure', width: 300, renderCell: (params : GridRenderCellParams) => params.value },
    { field: 'logP', headerName: 'LogP', width: 130 },
    { field: 'qed', headerName: 'QED Score', width: 130 },
    { field: 'mol_weight', headerName: 'Molecular Weight', width: 180 },
    { field: 'num_h_donors', headerName: 'H Donors', width: 130 },
  ];

  return (
    <>
      <div style={{ height: 400, width: '100%' }}>
      <DataGrid rows={rows} 
      columns={columns} 
      rowHeight={200} 
      checkboxSelection
      slots={{toolbar: GridToolbar }}/>
    </div>
    </>
  );
};

export default StructureOutput;
