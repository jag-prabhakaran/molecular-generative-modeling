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
  const [qedMax, setQedMax] = useState('0')
  const [qedMin, setQedMin] = useState('0')
  const [NumOfMolecules, setNumOfMolecules] = useState('5');
  const [selectedModel, setSelectedModel] = useState('scaffold-constrained');
  const [rationale, setRationale] = useState("OCc1cc[c:1]c(-c2ncccn2)c1")
  const [description, setModelDescription] = useState('')
  const [apiResponse, setAPIResponse] = useState(null);
  const [loading, setLoading] = useState(false);
  

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

  const handleQedMaxChange = (newQedMax) => {
    setQedMax(newQedMax);
  };

  const handleQedMinChange = (newQedMin) => {
    setQedMin(newQedMin)
  };

  const handleModelChange = (newModel) => {
    setSelectedModel(newModel);
  };
  


  const handleGenerate = async () => {
    var rationaleFixed = "OCc1cc[c:1]c(-c2ncccn2)c1";
    try {
      const data = {
        'model_type': selectedModel,
        'payload': {
          'scaffold_smile': smiles,
          'rationale': [rationaleFixed],
          'log_p_max': parseFloat(logPMin),
          'log_p_min': parseFloat(logPMax),
          'num_molecules': parseInt(NumOfMolecules),
          'qed_max': parseFloat(qedMax),
          'qed_min': parseFloat(qedMin)
        }
      };
      console.log('Sending data:', JSON.stringify(data));
      setLoading(true);
      const response = await (await fetch('https://ezu74lbfo2imcxnmbblg3hkhqq0oqtxo.lambda-url.us-east-1.on.aws/', {
        method: 'POST',
        body: JSON.stringify(data),
      })).json();
      setLoading(false);
      if(selectedModel === "multiobj-rationale"){
        setAPIResponse(response.filtered_output_objects);
      } else {
        setAPIResponse(response.filtered_smiles);
      }


      console.log('Reponse recieved successfully:', apiResponse);
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
        loading={loading}
        />
        <JSMEContainer 
        onSmilesChange={handleSmilesChange}
        />
        <Properties
        logPMax={logPMax}
        logPMin={logPMin}
        NumOfMolecules={NumOfMolecules}
        qedMax={qedMax}
        qedMin={qedMin}
        onLogPMaxChange={handleLogPMaxChange}
        onLogPMinChange={handleLogPMinChange}
        onNumOfMoleculesChange={handleNumOfMoleculesChange}
        onQedMaxChange={handleQedMaxChange}
        onQedMinChange={handleQedMinChange}
        />
        <div>
        <GenerateButton onGenerate={handleGenerate}/>
        {}
        </div>
        {apiResponse && <VisualOutput 
        apiResponse={apiResponse}
        />}
      </header>
    </div>
  );
};
  
export default App;