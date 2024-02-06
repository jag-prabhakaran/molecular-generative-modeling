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

  const MolRender = (props : { key: string, molecule: object})  => {
    const{key, molecule} = props
    return (
        <Box>
            molecule
        </Box>
    )
  }

  export default MolRender