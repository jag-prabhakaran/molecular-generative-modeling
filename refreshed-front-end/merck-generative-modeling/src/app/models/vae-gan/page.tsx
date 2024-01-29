"use client";
import React from 'react';
import MerckNavbar from '../../_components/MerckNavbar';
import { Box, Button, Card, CssBaseline, Toolbar, Typography } from '@mui/material';
import KetcherComponent from '@/app/_components/KetcherComponent';
import PropertyControls from '@/app/_components/PropertyControls';

const vaeGan: React.FC = () => {
    return (
        <Box className='flex'>
            <CssBaseline />
            <MerckNavbar />
            <Box component="main" className='flex flex-col flex-grow'>
                <Toolbar />
                <Box className='flex flex-row justify-evenly'>
                    <Box className='w-8/12 p-3' style={{height: "50vh"}}>
                        <KetcherComponent  />
                        {/* <Box className='w-full h-full' bgcolor="red"></Box> */}
                    </Box>
                    <Card className='w-3/12 p-3 overflow-scroll'>
                       <PropertyControls /> 
                    </Card>
                </Box>
                <Box className='flex flex-row justify-center'>
                    <Button variant="contained" >Generate Molecules</Button>
                </Box>
            </Box>
        </Box>
    );
};

export default vaeGan;
