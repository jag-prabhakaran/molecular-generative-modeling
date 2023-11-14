import React from 'react';
import ModelCheckpoint from './components/ModelCheckpoint';
import ModelDescription from './components/ModelDescription';
import JSMEContainer from './components/JSMEContainer';
import Properties from './components/Properties';
import GenerateButton from './components/GenerateButton';
import './App.css';


function App() {
  return (
    <div className="App">
      <header className="App-header">
        <ModelCheckpoint />
        <ModelDescription />
        <JSMEContainer />
        <Properties />
        <GenerateButton />
      </header>
    </div>
  );
}

export default App;
