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

interface PropertyControlProps {
    propertyValues: { [key: string]: string };
    handlePropertyChange: (property: string, value: string) => void;
}

const PropertyControls: React.FC<PropertyControlProps> = (props) => {
  const propertiesArray = [
    "logP Min",
    "logP Max",
    "num molecules",
    "qed Min",
    "qed Max",
    "more fake properties",
    "more fake properties",
    "more fake properties",
  ];


  const singlePropertyControl = (property: string) => (
    <Box className="flex flex-grow flex-row align-middle justify-center p-2">
      <TextField
        id={property}
        label={property}
        variant="outlined"
        value={props.propertyValues[property] || ""}
        onChange={(e) => props.handlePropertyChange(property, e.target.value)}
      />
      <Button size="small">Reset</Button>
    </Box>
  );

  return (
    <>{propertiesArray.map((property) => singlePropertyControl(property))}</>
  );
};

export default PropertyControls;
