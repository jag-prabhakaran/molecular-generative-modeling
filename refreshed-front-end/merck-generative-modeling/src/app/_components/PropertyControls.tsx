import { Box, Button, Card, CardActions, CardContent, Input, List, Slider, TextField, Typography } from '@mui/material';
import React from 'react';

interface PropertyControlProps{
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
        <Box className='flex flex-grow flex-row align-middle justify-center p-2'>
            {/* <CardContent> */}
                <TextField id="outlined-basic" label={property} variant="outlined" />
            {/* </CardContent> */}
            {/* <CardActions> */}
                <Button size="small">Reset</Button>
            {/* </CardActions> */}
        </Box>
        )
    return (
        <>
            {propertiesArray.map((property) => singlePropertyControl(property))}
        </>
    )
};

export default PropertyControls;
