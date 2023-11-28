import React, { useEffect, useState } from 'react';
import ModelCheckpoint from './components/ModelCheckpoint';
import ModelDescription from './components/ModelDescription';
import Properties from './components/Properties';
import GenerateButton from './components/GenerateButton';
import VisualOutput from './components/VisualOutput';
import './App.css';
import JSMEContainer from './components/JSMEContainer.js';


function App() {
  const [smiles, setSmiles] = useState("CC(C)(C(=O)O)c1ccc(cc1)C(O)CCCN2CCC(CC2)C(O)(*)c4ccccc4");
  const [logPMax, setLogPMax] = useState('0');
  const [logPMin, setLogPMin] = useState('0');
  const [NumOfMolecules, setNumOfMolecules] = useState('5');
  const [selectedModel, setSelectedModel] = useState('scaffold-constrained');
  const [rationale, setRationale] = useState("OCc1cc[c:1]c(-c2ncccn2)c1")
  const [description, setModelDescription] = useState('')
  const [VisualData, setVisualOuput] = useState(false)
  const [SMILES_LIST, setSMILES_LIST] = useState([])
  const [subSMILES_LIST, setsubSMILES_LIST] = useState([])
  

  const handleSmilesChange = (newSmiles) => {
    setSmiles(newSmiles);
    setRationale(newSmiles);
  };

  const handleLogPMaxChange = (newLogPMax) => {
    setLogPMax(newLogPMax);
  };

  const handleLogPMinChange = (newLogPMin) => {
    setLogPMin(newLogPMin)
  };

  const handleNumOfMoleculesChange = (newNumOfMolecules) => {
    setNumOfMolecules(newNumOfMolecules)
  };

  const handleModelChange = (newModel) => {
    setSelectedModel(newModel);
  };
  


  const handleGenerate = async () => {
    try {
      const data = {
        'model_type': selectedModel,
        'payload': {
          'scaffold_smile': smiles,
          'rationale': [rationale],
          'log_p_max': parseFloat(logPMin),
          'log_p_min': parseFloat(logPMax),
          'num_molecules': parseInt(NumOfMolecules)
        }
      };
      setSMILES_LIST([])
      setsubSMILES_LIST([])
      setVisualOuput(true)
      console.log('Data sent', data)
      const response = await fetch('https://ezu74lbfo2imcxnmbblg3hkhqq0oqtxo.lambda-url.us-east-1.on.aws/', {
        method: 'POST',
        body: JSON.stringify(data),
      });
  
      const result = await response.json();
      console.log('Data sent successfully:', result);
      
      switch(selectedModel) {
        case "scaffold-constrained":
          const smilesListScaffold = result.filtered_smiles.map(element => element.smile)
          const subsmileslistScaffold = result.filtered_smiles.map(element => "")
          setSMILES_LIST(smilesListScaffold)
          setsubSMILES_LIST(subsmileslistScaffold)
          break
        case "vae-gan":
          const smilesListGan = result.filtered_smiles.map(element => element.smile)
          const subsmileslistGan = result.filtered_smiles.map(element => "")
          setSMILES_LIST(smilesListGan)
          setsubSMILES_LIST(subsmileslistGan)
          break
        case "multiobj-rationale":
          const smileListRationale = result.filtered_output_objects.map(element => element.output_smile)
          const subsmilesRationale = result.filtered_output_objects.map(element => "")
          setSMILES_LIST(smileListRationale)
          setsubSMILES_LIST(subsmilesRationale)
          break
        default:
          return null;
      }
    } catch (error) {
      console.error('Error sending data:', error);
    }
  };

  useEffect (() => {
    switch(selectedModel) {
      case "scaffold-constrained":
        setModelDescription(`
        This model employs a SMILES-based Recurrent Neural Network (RNN) generative model to 
        achieve scaffold-constrained generation.
        The model takes as inputs :
        <ul>
            <li>a molecular scaffold in the form of a SMILES string with open positions marked by * (asterisks)
                [there can be one or many open positions in a single scaffold SMILES string]</li>
            <li>the upper bound of the number of molecules to be attempted to be generated. </li>
        </ul>
        The model outputs SMILES strings for the generated molecules.
        `);
        break;
      
      case "vae-gan":
        setModelDescription(`
        It's a simple GAN (Generative Adversarial Network) architecture that trains on small molecules 
        (number of atoms < 9) and outputs new (small) molecules from noise.
        `);
        break;
      
      case "multiobj-rationale":
        setModelDescription(`
        The model works by decomposing a molecule into a graph where each node corresponds to a specific substructure 
        of the molecule. This graph is then passed through a VAE (Variational Autoencoder) model to generate new molecules. 
        The user can specify a property and a input molecule(or molecules) and the model generates molecules from the input 
        molecules which improves that specific property. 
        Pro: User can mention a wide range of properties, has some of the best benchmarks. 
        Con: Doesn't allow user to control where the changes are made in the input structure of the molecule, 
        can be pretty slow compared to other molecules (takes ~0.4 seconds per molecule on a GPU)
        `)
        break;
      default:
        return null;
    }
  }, [selectedModel])

    
  
  
  return (
    <div className="App">
      <header className="App-header">
        <ModelCheckpoint 
        onModelChange={handleModelChange}
        />
        <ModelDescription 
        description={description}
        />
        <JSMEContainer 
        onSmilesChange={handleSmilesChange}
        />
        <Properties
        logPMax={logPMax}
        logPMin={logPMin}
        NumOfMolecules={NumOfMolecules}
        onLogPMaxChange={handleLogPMaxChange}
        onLogPMinChange={handleLogPMinChange}
        onNumOfMoleculesChange={handleNumOfMoleculesChange}
        />
        <GenerateButton onGenerate={handleGenerate}/>
        {VisualData && <VisualOutput 
        VisualData={VisualData} 
        SMILES_LIST={SMILES_LIST}
        subSMILES_LIST={subSMILES_LIST}
        />}
      </header>
    </div>
  );
};
  
export default App;