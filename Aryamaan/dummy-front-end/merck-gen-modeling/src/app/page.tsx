"use client";
import { useState } from "react";

export default function Home() {
  const [numMolecules, setNumMolecules] = useState(1);
  const [rationale, setRationale] = useState("OCc1cc[c:1]c(-c2ncccn2)c1");
  const [scaffold, setScaffold] = useState(
    "CC(C)(C(=O)O)c1ccc(cc1)C(O)CCCN2CCC(CC2)C(O)(*)c4ccccc4"
  );
  const [minLogP, setMinLogP] = useState(0);
  const [maxLogP, setMaxLogP] = useState(0);
  const [model, setModel] = useState("scaffold-constrained");
  const [fetching, setFetching] = useState(false);
  const [apiResponse, setApiResponse] = useState("");

  const handleNumMoleculesChange = (
    event: React.ChangeEvent<HTMLInputElement>
  ) => {
    setNumMolecules(Number(event.currentTarget.value));
  };

  const handleRationaleChange = (
    event: React.ChangeEvent<HTMLInputElement>
  ) => {
    setRationale(event.currentTarget.value);
  };

  const handleScaffoldChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setScaffold(event.currentTarget.value);
  };

  const handleMinLogPChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setMinLogP(Number(event.currentTarget.value));
  };

  const handleMaxLogPChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setMaxLogP(Number(event.currentTarget.value));
  };

  const handleModelChange = (event: React.ChangeEvent<HTMLSelectElement>) => {
    setModel(event.currentTarget.value);
  };

  const handleButtonClick = async (
    event: React.MouseEvent<HTMLButtonElement>
  ) => {
    try {
      setFetching(true);
      const response = await fetch(
        "https://ezu74lbfo2imcxnmbblg3hkhqq0oqtxo.lambda-url.us-east-1.on.aws/",
        {
          method: "POST",
          body: JSON.stringify({
            model_type: model,
            payload: {
              num_molecules: numMolecules,
              rationale: [rationale],
              scaffold_smile: scaffold,
              log_p_min: minLogP,
              log_p_max: maxLogP,
            },
          }),
        }
      );
      setFetching(false);
      if (!response.ok) {
        throw new Error("Network response was not ok");
      }

      const data = await response.json();
      setApiResponse(data);
    } catch (error) {
      console.error("Error:", error);
      setApiResponse("Error occurred during API call");
    }
  };

  const renderModelOptions = () => {
    switch (model) {
      case "scaffold-constrained":
        return (
          <div className="flex flex-col w-1/2">
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="logPMin">Log P Min:</label>
              <input
                onChange={handleMinLogPChange}
                className="border w-1/12"
                type="number"
                id="logPMin"
                name="log_p_min"
                value={minLogP}
                min="0"
              />
            </div>
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="logPMax">Log P Max:</label>
              <input
                onChange={handleMaxLogPChange}
                className="border w-1/12"
                type="number"
                id="logPMax"
                name="log_p_max"
                value={maxLogP}
                min="1"
              />
            </div>
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="scaffoldSmile">Scaffold SMILE:</label>
              <input
                onChange={handleScaffoldChange}
                className="border w-full"
                type="text"
                id="scaffoldSmile"
                name="scaffold_smile"
                value={scaffold}
              />
            </div>
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="numMolecules">Number of Molecules:</label>
              <input
                onChange={handleNumMoleculesChange}
                className="border w-1/12"
                type="number"
                id="numMolecules"
                name="num_molecules"
                value={numMolecules}
                min="1"
              />
            </div>
          </div>
        );
      case "vae-gan":
        return (
          <div className="flex flex-col w-1/2">
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="logPMin">Log P Min:</label>
              <input
                onChange={handleMinLogPChange}
                className="border w-1/12"
                type="number"
                id="logPMin"
                name="log_p_min"
                value={minLogP}
                min="0"
              />
            </div>
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="logPMax">Log P Max:</label>
              <input
                onChange={handleMaxLogPChange}
                className="border w-1/12"
                type="number"
                id="logPMax"
                name="log_p_max"
                value={maxLogP}
                min="1"
              />
            </div>
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="numMolecules">Number of Molecules:</label>
              <input
                onChange={handleNumMoleculesChange}
                className="border w-1/12"
                type="number"
                id="numMolecules"
                name="num_molecules"
                value={numMolecules}
              />
            </div>
          </div>
        );
      case "multiobj-rationale":
        return (
          <div className="flex flex-col w-1/2">
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="logPMin">Log P Min:</label>
              <input
                onChange={handleMinLogPChange}
                className="border w-1/12"
                type="number"
                id="logPMin"
                name="log_p_min"
                value={minLogP}
                min="0"
              />
            </div>
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="logPMax">Log P Max:</label>
              <input
                onChange={handleMaxLogPChange}
                className="border w-1/12"
                type="number"
                id="logPMax"
                name="log_p_max"
                value={maxLogP}
                min="1"
              />
            </div>
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="rationaleString">Rationale String</label>
              <input
                onChange={handleRationaleChange}
                className="border w-full"
                type="text"
                id="rationaleString"
                name="scaffold_smile"
                value={rationale}
              />
            </div>
            <div className="flex flex-row gap-2 justify-between">
              <label htmlFor="numMolecules">Number of Molecules:</label>
              <input
                className="border w-1/12"
                type="number"
                id="numMolecules"
                name="num_molecules"
                value={numMolecules}
                min="1"
              />
            </div>
          </div>
        );
      default:
        return null;
    }
  };

  return (
    <div className="p-5">
      <div className="flex flex-col items-center gap-3">
        <h1 className="text-3xl font-bold">{model}</h1>
        <div className="flex flex-row items-center gap-3 w-full justify-evenly">
          <select className="border-2" onChange={handleModelChange}>
            <option value="scaffold-constrained">scaffold-constrained</option>
            <option value="vae-gan">vae-gan</option>
            <option value="multiobj-rationale">multiobj-rationale</option>
          </select>
          {renderModelOptions()}
          <button
            type="button"
            className="text-white 
        bg-blue-700 
        hover:bg-blue-800 
        focus:outline-none 
        focus:ring-4 
        focus:ring-blue-300 
        font-medium 
        rounded-full 
        text-sm 
        px-5 
        py-2.5 
        text-center 
        me-2 
        mb-2 
        dark:bg-blue-600 
        dark:hover:bg-blue-700 
        dark:focus:ring-blue-800
        my-3"
            onClick={handleButtonClick}
          >
            {fetching ? (
              <svg
                aria-hidden="true"
                role="status"
                className="inline w-6 h-6 me-3 text-white animate-spin"
                viewBox="0 0 100 101"
                fill="none"
                xmlns="http://www.w3.org/2000/svg"
              >
                <path
                  d="M100 50.5908C100 78.2051 77.6142 100.591 50 100.591C22.3858 100.591 0 78.2051 0 50.5908C0 22.9766 22.3858 0.59082 50 0.59082C77.6142 0.59082 100 22.9766 100 50.5908ZM9.08144 50.5908C9.08144 73.1895 27.4013 91.5094 50 91.5094C72.5987 91.5094 90.9186 73.1895 90.9186 50.5908C90.9186 27.9921 72.5987 9.67226 50 9.67226C27.4013 9.67226 9.08144 27.9921 9.08144 50.5908Z"
                  fill="#E5E7EB"
                />
                <path
                  d="M93.9676 39.0409C96.393 38.4038 97.8624 35.9116 97.0079 33.5539C95.2932 28.8227 92.871 24.3692 89.8167 20.348C85.8452 15.1192 80.8826 10.7238 75.2124 7.41289C69.5422 4.10194 63.2754 1.94025 56.7698 1.05124C51.7666 0.367541 46.6976 0.446843 41.7345 1.27873C39.2613 1.69328 37.813 4.19778 38.4501 6.62326C39.0873 9.04874 41.5694 10.4717 44.0505 10.1071C47.8511 9.54855 51.7191 9.52689 55.5402 10.0491C60.8642 10.7766 65.9928 12.5457 70.6331 15.2552C75.2735 17.9648 79.3347 21.5619 82.5849 25.841C84.9175 28.9121 86.7997 32.2913 88.1811 35.8758C89.083 38.2158 91.5421 39.6781 93.9676 39.0409Z"
                  fill="currentColor"
                />
              </svg>
            ) : null}
            Make Request
          </button>
        </div>
        <div className="flex flex-col">
          <strong>API Response:</strong>
          {fetching ? <p>Fetching...</p> : null}
          <p>{JSON.stringify(apiResponse)}</p>
        </div>
      </div>
    </div>
  );
}
