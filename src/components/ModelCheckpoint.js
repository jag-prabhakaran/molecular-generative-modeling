// ModelCheckpoint.js
import React from 'react';

const ModelCheckpoint = () => {
  return (
    <div>
        <div class="ModelCheckpointText" style={{
            color: 'black',
            fontSize: '3.5vh',
            fontFamily: 'monospace',
            fontWeight: 400,
            top: '5vh',
            left: '2vw',
            wordWrap: 'break-word',
            position: 'absolute'
            }}>Model Checkpoint
        </div>
        <div class="ModelSelection">
            <select style={{
                width: '243px',
                height: '43px',
                backgroundColor: '#D0C1C1',
                borderRadius: '10px',
                top: '11vh',
                left: '2vw',
                position: 'absolute'
            }}>
                <option value="ScaffoldConstraint">SAMOA</option>
                <option value="VAEGAN">Graph Translation</option>
                <option value="MolGan">MolGan</option>
            </select>
        </div>
    </div>
  );
};

export default ModelCheckpoint;
