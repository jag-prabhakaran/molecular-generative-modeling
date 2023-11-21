//import React, { useState } from 'react';
import ModelCheckpoint from './components/ModelCheckpoint';
import ModelDescription from './components/ModelDescription';
import JSMEContainer from './components/JSMEContainer';
import Properties from './components/Properties';
import GenerateButton from './components/GenerateButton';
import './App.css';


function App() {
  // const [smiles, setSmiles] = useState('');
  // const [logPMax, setLogPMax] = useState('0');
  // const [logPMin, setLogPMin] = useState('0');
  // const [NumOfMolecules, setNumOfMolecules] = useState('5');
  // const [selectedModel, setSelectedModel] = useState('scaffold-constrained');

  // const handleSmilesChange = (newSmiles) => {
  //   setSmiles(newSmiles);
  // };

  // const handleLogPMaxChange = (newLogPMax) => {
  //   setLogPMax(newLogPMax);
  // };

  // const handleLogPMinChange = (newLogPMin) => {
  //   setLogPMin(newLogPMin)
  // };

  // const handleNumOfMoleculesChange = (newNumOfMolecules) => {
  //   setNumOfMolecules(newNumOfMolecules)
  // };

  // const handleModelChange = (newModel) => {
  //   setSelectedModel(newModel);
  // };

  const handleGenerate = async () => {
    try {
      const data = {
        'model_type': "vae-gan",
        'payload': {
          'scaffold_smile': 'CC(C)(C(=O)O)c1ccc(cc1)C(O)CCCN2CCC(CC2)C(O)(*)c4ccccc4',
          'log_p_max': 10,
          'log_p_min': 1,
          'num_molecules': 10
        }
      };
  
      const response = await fetch('https://ezu74lbfo2imcxnmbblg3hkhqq0oqtxo.lambda-url.us-east-1.on.aws/', {
        method: 'POST',
        body: JSON.stringify(data),
      });
  
      const result = await response.json();
      console.log('Data sent successfully:', result);
    } catch (error) {
      console.error('Error sending data:', error);
    }
  };
  
  return (
    <div className="App">
      <header className="App-header">
        <ModelCheckpoint 
        //onModelChange={handleModelChange}
        />
        <ModelDescription />
        <JSMEContainer 
        //onSmilesChange={handleSmilesChange}
        />
        <Properties
        // logPMax={logPMax}
        // logPMin={logPMin}
        // NumOfMolecules={NumOfMolecules}
        // onLogPMaxChange={handleLogPMaxChange}
        // onLogPMinChange={handleLogPMinChange}
        // onNumOfMoleculesChange={handleNumOfMoleculesChange}
        />
        <GenerateButton onGenerate={handleGenerate}/>
      </header>
    </div>
  );
};
  
export default App;