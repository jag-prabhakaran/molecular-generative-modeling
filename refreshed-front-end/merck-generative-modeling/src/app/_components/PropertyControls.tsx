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

/*interface PropertyControlProps {
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
    "solubility @ ph7",
    "FaSSIF solubility",
    "num H donors",
    "HPLC logD",
    "membrane permeability",
    "blood/brain barrier permeability"
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
};*/

interface PropertyControlProps {
  propertyValues: { [key: string]: string };
  handlePropertyChange: (property: string, value: string) => void;
}

const PropertyControls = (props: PropertyControlProps) => {
  const textFieldProperties = ["num molecules"];

  const sliderProperties = ["logP Min", "logP Max", "qed Min", "qed Max"];

  const combinedProperties = [...textFieldProperties, ...sliderProperties];

  const singlePropertyControl = (property: string) => {
    // Check if the current property should use a slider
    if (sliderProperties.includes(property)) {
      return (
        <Box className="flex flex-grow flex-row align-middle justify-center p-2">
          <Typography id={`slider-${property}`} gutterBottom>
            {property}
          </Typography>
          <Box className="w-1/2 mx-2">
            <Slider
              value={
                props.propertyValues[property]
                  ? Number(props.propertyValues[property])
                  : 0
              }
              onChange={(e, value) =>
                props.handlePropertyChange(property, String(value))
              }
              aria-labelledby={`slider-${property}`}
              step={0.1}
              min={0}
              max={property == "qed Min" || property == "qed Max" ? 1 : 10}
              valueLabelDisplay="auto"
            />
          </Box>
          <Input
            value={
              props.propertyValues[property]
                ? Number(props.propertyValues[property])
                : 5
            }
            size="small"
            onChange={(e) =>
              props.handlePropertyChange(property, e.target.value)
            }
            inputProps={{
              step: 1,
              min: 0,
              max: 10,
              type: "number",
              "aria-labelledby": "slider-" + property,
            }}
          />
        </Box>
      );
    } else {
      // If not a slider property, return a TextField
      return (
        <Box className="flex flex-grow flex-row align-middle justify-center p-2">
          <TextField
            id={property}
            label={property}
            variant="outlined"
            value={props.propertyValues[property] || ""}
            onChange={(e) =>
              props.handlePropertyChange(property, e.target.value)
            }
          />
          <Button size="small">Reset</Button>
        </Box>
      );
    }
  };

  return (
    <>{combinedProperties.map((property) => singlePropertyControl(property))}</>
  );
};

export default PropertyControls;
