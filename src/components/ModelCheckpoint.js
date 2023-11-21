// ModelCheckpoint.js

// ModelCheckpoint.js
import React from 'react';

const ModelCheckpoint = ({ onModelChange }) => {
  const handleModelChange = (event) => {
    const selectedModel = event.target.value;
    onModelChange(selectedModel);
  };

  return (
    <div>
      <div className="ModelCheckpointText" style={{
        color: 'black',
        fontSize: '3.5vh',
        fontFamily: 'monospace',
        fontWeight: 400,
        top: '5vh',
        left: '2vw',
        wordWrap: 'break-word',
        position: 'absolute',
      }}>
        Model Checkpoint
      </div>
      <div className="ModelSelection">
        <select
          onChange={handleModelChange}
          style={{
            width: '243px',
            height: '43px',
            backgroundColor: '#D0C1C1',
            borderRadius: '10px',
            top: '11vh',
            left: '2vw',
            position: 'absolute',
          }}
        >
          <option value="scaffold-constrained">SAMOA</option>
          <option value="vae-gan">Graph Translation</option>
          <option value="MolGan">MolGan</option>
        </select>
      </div>
    </div>
  );
};

export default ModelCheckpoint;