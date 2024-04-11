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
import {
  DataGrid,
  GridColDef,
  GridRenderCellParams,
  GridToolbar,
} from "@mui/x-data-grid";
import MolRender from "./MolRender";
import MoleculeStructure from "./MoleculeStructure";

const StructureOutput: React.FC<{
  response: any;
  isMultiObj: boolean;
  input_smile: string;
}> = ({ response, isMultiObj, input_smile }) => {
  const round = (value: number) => Math.round(value * 1000) / 1000;

  const rows = Object.keys(response).map((key) => ({
    id: key,
    smile: key,
    logP: round(response[key].logP),
    qed: round(response[key].qed),
    mol_weight: round(response[key].mol_weight),
    num_h_donors: round(response[key].num_h_donors),
    moleculeStructure: (
      <MoleculeStructure
        structure={key}
        id={`mol_${key}`}
        subStructure={input_smile}
        svgMode={true}
      ></MoleculeStructure>
    ),
  }));

  const columns = [
    {
      field: "moleculeStructure",
      headerName: "Structure",
      width: 300,
      renderCell: (params: GridRenderCellParams) => params.value,
      disableExport: true
    },
    { field: "smile", headerName: "SMILE", width: 130, hide: true},
    { field: "logP", headerName: "LogP", width: 130 },
    { field: "qed", headerName: "QED Score", width: 130 },
    { field: "mol_weight", headerName: "Molecular Weight", width: 180 },
    { field: "num_h_donors", headerName: "H Donors", width: 130 },
  ];

  const columnVisibilityModel = {smile:false}
  return (
    <>
      <div style={{ height: 400, width: "100%" }}>
        <DataGrid
          rows={rows}
          columns={columns}
          rowHeight={200}
          checkboxSelection
          slots={{ toolbar: GridToolbar }}
          slotProps={{ toolbar: {printOptions: {disableToolbarButton: true},
          csvOptions: {allColumns: true} }}}
          columnVisibilityModel={columnVisibilityModel}
        />
      </div>
    </>
  );
};

export default StructureOutput;
