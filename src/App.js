import React from 'react';
import logo from './logo.svg';
import MyComponent from './components/MyComponent'; // Make sure to adjust the path if needed
import './App.css';

function App() {
  return (
    <div className="App">
      <header className="App-header">
        <MyComponent /> {/* Add your MyComponent here */}
      </header>
    </div>
  );
}

export default App;
