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
  propertyValues: { [key: string]: string }
  handlePropertyChange: (property: string, value: string) => void;
}

const PropertyControls = (props: PropertyControlProps) => {
  const textFieldProperties = ["upper bound"];

  const sliderProperties = ["logP", "qed"];

  const combinedProperties = [...textFieldProperties, ...sliderProperties];

  const singlePropertyControl = (property: string) => {
    // Check if the current property should use a slider
    if (sliderProperties.includes(property)) {
      return (
        <Box className="flex flex-grow flex-col align-middle justify-center p-2">
          <Typography id={`slider-${property}`} gutterBottom>
            {property}
          </Typography>
          <Box className="flex flex-row align middle justify-center p-2">
          <Input
            className="pr-5"
            value={props.propertyValues[`${property} Min`] || 0}
            size="small"
            onChange={(e) =>
              props.handlePropertyChange(`${property} Min`, e.target.value)
            }
            inputProps={{
              step: 0.1,
              min: 0,
              max: property === "qed" ? 1 : 10,
              type: "number",
              "aria-labelledby": "input-" + property + "min",
            }}
          />
          <Box className="w-1/2 mx-2">
            <Slider
              value={[
                props.propertyValues[`${property} Min`]
                  ? Number(props.propertyValues[`${property} Min`])
                  : 0,
                props.propertyValues[`${property} Max`]
                  ? Number(props.propertyValues[`${property} Max`])
                  : (property === 'qed' ? 1 : 10)
              ]}
              onChange={(e, value) => {
                const [min,max] = value as number[]
                props.handlePropertyChange(`${property} Min`, String(min))
                props.handlePropertyChange(`${property} Max`, String(max))
              }}
              aria-labelledby={`slider-${property}`}
              step={0.1}
              min={0}
              max={property === "qed" ? 1 : 10}
              valueLabelDisplay="auto"
              disableSwap
            />
          </Box>
          <Input
            className="pl-5"
            value={props.propertyValues[`${property} Max`] || (property === 'qed' ? 1 : 10)}
            size="small"
            onChange={(e) =>
              props.handlePropertyChange(`${property} Max`, e.target.value)
            }
            inputProps={{
              step: 0.1,
              min: 0,
              max: property === "qed" ? 1 : 10,
              type: "number",
              "aria-labelledby": "input-" + property + "max"
            }}
          />
          </Box>
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
