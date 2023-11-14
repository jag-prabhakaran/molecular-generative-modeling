import React, { useState } from 'react';
import ModelCheckpoint from './components/ModelCheckpoint';
import ModelDescription from './components/ModelDescription';
import JSMEContainer from './components/JSMEContainer';
import Properties from './components/Properties';
import GenerateButton from './components/GenerateButton';
import VisualOutput from './components/VisualOutput'
import './App.css';


function App() {
  const [smiles, setSmiles] = useState('');
  const [logP, setLogP] = useState('');
  const [NumOfMolecules, setNumOfMolecules] = useState('');
  const [visualData, setVisualData] = useState(false);

  const handleSmilesChange = (newSmiles) => {
    setSmiles(newSmiles);
  };

  const handleLogPChange = (newLogP) => {
    setLogP(newLogP);
  };

  const handleNumOfMoleculesChange = (newNumOfMolecules) => {
    setNumOfMolecules(newNumOfMolecules)
  };

  const handleGenerate = async () => {
    const data = {
      smiles: smiles,
      logP: logP,
      NumOfMolecules: NumOfMolecules
    };
    setVisualData(true)
    console.log('Data to be sent:', data)
  }
  return (
    <div className="App">
      <header className="App-header">
        <ModelCheckpoint />
        <ModelDescription />
        <JSMEContainer 
        onSmilesChange={handleSmilesChange}
        />
        <Properties
        logP={logP}
        NumOfMolecules={NumOfMolecules}
        onLogPChange={handleLogPChange}
        onNumOfMoleculesChange={handleNumOfMoleculesChange}
        />
        <GenerateButton onGenerate={handleGenerate}/>
        {visualData && < VisualOutput />}
      </header>
    </div>
  );
}

export default App;
