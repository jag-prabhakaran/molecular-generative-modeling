import React from 'react';
import ModelCheckpoint from './components/ModelCheckpoint';
import ModelDescription from './components/ModelDescription';
import JSMEContainer from './components/JSMEContainer';
import './App.css';


function App() {
  return (
    <div className="App">
      <header className="App-header">
        <ModelCheckpoint />
        <ModelDescription />
        <JSMEContainer />
      </header>
    </div>
  );
}

export default App;
