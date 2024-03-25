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
    smile: molecule.smile,
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
    { field: 'moleculeStructure', 
    headerName: 'Structure', 
    width: 300,
    disableExport: true, 
    renderCell: (params : GridRenderCellParams) => params.value },

    {field : 'smile',
    headerName: 'SMILE',
    width: 300,
    hide: true},

    { field: 'logP', 
    headerName: 'LogP',
    width: 130,
    type : 'number' as const },

    { field: 'qed',
     headerName: 'QED Score',
      width: 130, 
     type: 'number' as const},

    { field: 'mol_weight',
     headerName: 'Molecular Weight',
      width: 180, 
     type: 'number' as const},

    { field: 'num_h_donors',
     headerName: 'H Donors',
      width: 130 ,
    type: 'number' as const},
  ];

  const columnVisibilityModel = {smile: false}
  return (
    <>
      <div style={{ height: 400, width: '100%' }}>
      <DataGrid rows={rows} 
      columns={columns} 
      rowHeight={200} 
      checkboxSelection
      slots={{toolbar: GridToolbar }}
      slotProps={{ toolbar: { printOptions: { disableToolbarButton: true },
                              csvOptions: {allColumns: true}} }}
      columnVisibilityModel={columnVisibilityModel}
      />
    </div>
    </>
  );
};

export default StructureOutput;
