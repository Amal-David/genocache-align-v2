"""Simple seeding for Chr22 index"""
import torch
import numpy as np
import faiss

class Chr22Seeding:
    def __init__(self, model_path, index_path, positions_path, seed_len=512, K=32, device='cuda'):
        self.seed_len = seed_len
        self.K = K
        self.device = device
        
        # Load model
        print(f"Loading model from {model_path}...")
        from encoder_nal import NALEncoder
        checkpoint = torch.load(model_path, map_location=device)
        self.model = NALEncoder(
            emb_dim=checkpoint.get('emb_dim', 128),
            seed_len=checkpoint.get('seed_len', seed_len),
            vocab_size=5,
            hidden_dim=128,
            num_layers=4
        ).to(device)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model.eval()
        print(f"  ✅ Loaded (batch {checkpoint['batch']}, loss {checkpoint['loss']:.4f})")
        
        # Load index
        print(f"Loading index from {index_path}...")
        self.index = faiss.read_index(index_path)
        print(f"  ✅ {self.index.ntotal:,} vectors")
        
        # Load positions
        print(f"Loading positions from {positions_path}...")
        pos_data = np.load(positions_path, allow_pickle=True)
        self.chrs = pos_data['chr']
        self.positions = pos_data['pos']
        self.strands = pos_data['strand']
        print(f"  ✅ {len(self.positions):,} positions")
        
        self.base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    
    def get_anchors(self, read_seq, num_seeds=6, K=32):
        """Get anchors for a read"""
        anchors = []
        read_len = len(read_seq)
        
        # Extract seeds uniformly
        for i in range(num_seeds):
            start_pos = int(i * (read_len - self.seed_len) / max(1, num_seeds - 1)) if num_seeds > 1 else 0
            seed_seq = read_seq[start_pos:start_pos + self.seed_len]
            
            # Convert to indices
            seed_idx = [self.base_to_idx.get(b, 4) for b in seed_seq]
            seed_tensor = torch.tensor([seed_idx], dtype=torch.long).to(self.device)
            
            # Encode
            with torch.no_grad():
                embedding = self.model(seed_tensor).cpu().numpy()
            
            # Search index
            faiss.normalize_L2(embedding)
            similarities, indices = self.index.search(embedding, K)
            
            # Create anchors
            for j, (sim, idx) in enumerate(zip(similarities[0], indices[0])):
                anchors.append({
                    'seed_idx': i,
                    'seed_pos': start_pos,  # Position in read
                    'ref_chr': str(self.chrs[idx]),
                    'ref_pos': int(self.positions[idx]),
                    'strand': str(self.strands[idx]),
                    'similarity': float(sim),
                    'rank': j
                })
        
        return anchors
