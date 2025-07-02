#!/usr/bin/env python3
"""
Corrected RD Scanner for mosdepth output - with improved coverage calculation and filtering
"""

import argparse
import numpy as np
import pandas as pd
import gzip
from pathlib import Path
import logging
from typing import List, Tuple, Optional
import time


class CorrectedMosdepthRDScanner:
    def __init__(self, coverage_file: str, reference_genome_size: int = 4411532):
        """Initialize Corrected RD Scanner"""
        self.coverage_file = coverage_file
        self.reference_genome_size = reference_genome_size
        self.logger = logging.getLogger(__name__)
        
        # Default parameters
        self.coverage_threshold_factor = 0.1
        self.merge_distance = 1500
        self.min_rd_size = 500
        self.max_rd_size = 30000
        self.known_rd_threshold = 0.05
        
    def load_mosdepth_coverage(self) -> pd.DataFrame:
        """Load coverage data efficiently"""
        self.logger.info(f"Loading coverage data from {self.coverage_file}")
        
        if self.coverage_file.endswith('.gz'):
            df = pd.read_csv(self.coverage_file, sep='\t', compression='gzip',
                           names=['chrom', 'start', 'end', 'coverage'])
        else:
            df = pd.read_csv(self.coverage_file, sep='\t',
                           names=['chrom', 'start', 'end', 'coverage'])
        
        # Validate data
        df['bases'] = df['end'] - df['start']
        self.logger.info(f"Loaded {len(df)} coverage records, covering {df['bases'].sum()} bases")
        
        return df
    
    def find_low_coverage_regions(self, coverage_df: pd.DataFrame) -> List[Tuple[str, int, int]]:
        """Find regions with low coverage using improved filtering"""
        # Calculate average coverage properly
        total_bases = coverage_df['bases'].sum()
        weighted_cov = (coverage_df['bases'] * coverage_df['coverage']).sum()
        avg_coverage = weighted_cov / total_bases
        threshold = avg_coverage * self.coverage_threshold_factor
        
        self.logger.info(f"Average coverage: {avg_coverage:.2f}")
        self.logger.info(f"Coverage threshold: {threshold:.2f}")
        
        # Find truly low coverage regions
        low_coverage_df = coverage_df[coverage_df['coverage'] < threshold].copy()
        
        if low_coverage_df.empty:
            return []
        
        # Merge nearby regions and filter by actual coverage
        regions = []
        for chrom, group in low_coverage_df.groupby('chrom'):
            group = group.sort_values('start')
            merged = self._merge_regions_validated(group, coverage_df, chrom, threshold)
            regions.extend([(chrom, start, end) for start, end in merged])
        
        return regions
    
    def _merge_regions_validated(self, group: pd.DataFrame, full_coverage_df: pd.DataFrame, 
                                chrom: str, threshold: float) -> List[Tuple[int, int]]:
        """Merge regions with validation of actual coverage"""
        if group.empty:
            return []
        
        merged = []
        current_start = None
        current_end = None
        
        for _, row in group.iterrows():
            if current_start is None:
                current_start = row['start']
                current_end = row['end']
            elif row['start'] - current_end <= self.merge_distance:
                current_end = max(current_end, row['end'])
            else:
                # Validate the merged region before adding
                if current_end - current_start >= self.min_rd_size:
                    # Check actual coverage in the merged region
                    region_cov = self._calculate_exact_coverage(
                        full_coverage_df, chrom, current_start, current_end
                    )
                    if region_cov < threshold:  # Only add if truly low coverage
                        merged.append((current_start, current_end))
                
                current_start = row['start']
                current_end = row['end']
        
        # Don't forget the last region
        if current_start is not None and current_end - current_start >= self.min_rd_size:
            region_cov = self._calculate_exact_coverage(
                full_coverage_df, chrom, current_start, current_end
            )
            if region_cov < threshold:
                merged.append((current_start, current_end))
        
        return merged
    
    def _calculate_exact_coverage(self, coverage_df: pd.DataFrame, 
                                 chrom: str, start: int, end: int) -> float:
        """Calculate exact coverage for a region"""
        # Get all segments that overlap with our region
        overlapping = coverage_df[
            (coverage_df['chrom'] == chrom) & 
            (coverage_df['end'] > start) & 
            (coverage_df['start'] < end)
        ].copy()
        
        if overlapping.empty:
            return 0.0
        
        # Adjust boundaries for partial overlaps
        overlapping['adjusted_start'] = overlapping['start'].clip(lower=start)
        overlapping['adjusted_end'] = overlapping['end'].clip(upper=end)
        overlapping['adjusted_bases'] = overlapping['adjusted_end'] - overlapping['adjusted_start']
        
        # Calculate weighted average
        total_bases = overlapping['adjusted_bases'].sum()
        weighted_cov = (overlapping['adjusted_bases'] * overlapping['coverage']).sum()
        
        return weighted_cov / total_bases if total_bases > 0 else 0.0
    
    def analyze_region_coverage(self, coverage_df: pd.DataFrame, 
                              chrom: str, start: int, end: int, 
                              chrom_avg_cache: Optional[dict] = None) -> dict:
        """Analyze coverage for a single region with caching"""
        # Use exact coverage calculation
        avg_coverage = self._calculate_exact_coverage(coverage_df, chrom, start, end)
        
        # Calculate or use cached chromosome average
        if chrom_avg_cache and chrom in chrom_avg_cache:
            chrom_coverage = chrom_avg_cache[chrom]
        else:
            chrom_data = coverage_df[coverage_df['chrom'] == chrom]
            total_bases = chrom_data['bases'].sum()
            weighted_cov = (chrom_data['bases'] * chrom_data['coverage']).sum()
            chrom_coverage = weighted_cov / total_bases if total_bases > 0 else 1
            
            if chrom_avg_cache is not None:
                chrom_avg_cache[chrom] = chrom_coverage
        
        proportion = avg_coverage / chrom_coverage if chrom_coverage > 0 else 0
        
        return {
            'coverage': avg_coverage,
            'chrom_coverage': chrom_coverage,
            'proportion': proportion
        }
    
    def categorize_deletion(self, proportion: float) -> str:
        """Categorize deletion based on coverage proportion"""
        if proportion < 0.01:
            return "Complete deletion"
        elif proportion < 0.05:
            return "Strong deletion"
        elif proportion < 0.10:
            return "Moderate deletion"
        else:
            return "Weak signal"
    
    def analyze_novel_rds(self, coverage_df: pd.DataFrame) -> pd.DataFrame:
        """Find and analyze novel RD regions with improved validation"""
        self.logger.info("Detecting novel RD regions...")
        
        # Find low coverage regions
        low_coverage_regions = self.find_low_coverage_regions(coverage_df)
        
        # Filter by size
        filtered_regions = [(c, s, e) for c, s, e in low_coverage_regions 
                          if self.min_rd_size <= (e - s) <= self.max_rd_size]
        
        self.logger.info(f"Found {len(filtered_regions)} candidate regions after size filtering")
        
        # Analyze each region
        results = []
        chrom_avg_cache = {}
        
        for chrom, start, end in filtered_regions:
            analysis = self.analyze_region_coverage(coverage_df, chrom, start, end, chrom_avg_cache)
            
            # Additional validation - skip regions with high coverage
            if analysis['proportion'] >= self.coverage_threshold_factor:
                self.logger.warning(f"Skipping region {chrom}:{start}-{end} with proportion {analysis['proportion']:.3f} >= threshold")
                continue
            
            results.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'size': end - start,
                'type': 'DEL',
                'coverage': analysis['coverage'],
                'chrom_coverage': analysis['chrom_coverage'],
                'proportion': analysis['proportion'],
                'deletion_category': self.categorize_deletion(analysis['proportion']),
                'known_rd': ''
            })
        
        self.logger.info(f"Final: {len(results)} validated novel RD regions")
        return pd.DataFrame(results)
    
    def analyze_known_rds(self, coverage_df: pd.DataFrame, known_rds_file: str) -> pd.DataFrame:
        """Analyze coverage in known RD regions with improved accuracy"""
        self.logger.info(f"Analyzing known RDs from {known_rds_file}")
        
        # Load known RDs
        known_rds = pd.read_csv(known_rds_file, sep='\t', 
                              names=['chrom', 'start', 'end', 'name'])
        
        # Remove H37Rv reference entry
        known_rds = known_rds[known_rds['name'] != 'H37Rv'].reset_index(drop=True)
        
        results = []
        chrom_avg_cache = {}
        
        for _, rd in known_rds.iterrows():
            analysis = self.analyze_region_coverage(
                coverage_df, rd['chrom'], rd['start'], rd['end'], chrom_avg_cache
            )
            
            # Log significant differences for debugging
            if rd['name'] == 'RD207':
                self.logger.info(f"RD207 analysis: coverage={analysis['coverage']:.2f}, proportion={analysis['proportion']:.3f}")
            
            results.append({
                'RD': rd['name'],
                'Chromosome': rd['chrom'],
                'Start': rd['start'],
                'End': rd['end'],
                'Size': rd['end'] - rd['start'],
                'Coverage': analysis['coverage'],
                'ChromCoverage': analysis['chrom_coverage'],
                'Proportion': analysis['proportion'],
                'Status': 'Present' if analysis['proportion'] > self.known_rd_threshold else 'Absent'
            })
        
        return pd.DataFrame(results)
    
    def annotate_with_known_rds(self, novel_rds: pd.DataFrame, known_rds_file: str) -> pd.DataFrame:
        """Annotate novel RDs with known RD overlaps"""
        if novel_rds.empty or not Path(known_rds_file).exists():
            return novel_rds
        
        known_rds = pd.read_csv(known_rds_file, sep='\t', 
                              names=['chrom', 'start', 'end', 'name'])
        
        for idx, row in novel_rds.iterrows():
            overlapping = []
            
            # Find overlapping known RDs
            chrom_rds = known_rds[known_rds['chrom'] == row['chrom']]
            
            for _, rd in chrom_rds.iterrows():
                # Check overlap
                if row['end'] > rd['start'] and row['start'] < rd['end']:
                    # Check size similarity
                    known_size = rd['end'] - rd['start']
                    novel_size = row['size']
                    ratio = novel_size / known_size
                    
                    if 0.8 <= ratio <= 1.4:
                        overlapping.append(rd['name'])
            
            novel_rds.at[idx, 'known_rd'] = '|'.join(overlapping)
        
        return novel_rds
    
    def run(self, known_rds_file: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Run the complete RD detection pipeline"""
        start_time = time.time()
        
        # Load coverage data
        coverage_df = self.load_mosdepth_coverage()
        self.logger.info(f"Loaded {len(coverage_df)} coverage records in {time.time() - start_time:.2f}s")
        
        # Find novel RDs
        novel_rds = self.analyze_novel_rds(coverage_df)
        
        # Annotate with known RDs if provided
        if known_rds_file:
            novel_rds = self.annotate_with_known_rds(novel_rds, known_rds_file)
        
        # Analyze known RDs
        known_rds_analysis = pd.DataFrame()
        if known_rds_file:
            known_rds_analysis = self.analyze_known_rds(coverage_df, known_rds_file)
        
        self.logger.info(f"Analysis completed in {time.time() - start_time:.2f}s")
        
        return novel_rds, known_rds_analysis


def main():
    parser = argparse.ArgumentParser(description='Corrected RD detection from mosdepth coverage')
    parser.add_argument('coverage_file', help='Mosdepth per-base coverage file')
    parser.add_argument('-k', '--known-rds', required=True, help='BED file with known RD regions')
    parser.add_argument('-n', '--output-novel', default='novel_rds.tsv',
                       help='Output file for novel RDs')
    parser.add_argument('-o', '--output-known', default='known_rds_analysis.tsv',
                       help='Output file for known RDs analysis')
    parser.add_argument('-c', '--coverage-threshold', type=float, default=0.1,
                       help='Coverage threshold factor (default: 0.1)')
    parser.add_argument('-r', '--known-threshold', type=float, default=0.05,
                       help='Threshold for known RD presence (default: 0.05)')
    parser.add_argument('--min-size', type=int, default=500,
                       help='Minimum RD size (default: 500)')
    parser.add_argument('--max-size', type=int, default=30000,
                       help='Maximum RD size (default: 30000)')
    parser.add_argument('--merge-distance', type=int, default=1500,
                       help='Distance to merge nearby regions (default: 1500)')
    parser.add_argument('--log-level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level (default: INFO)')
    
    args = parser.parse_args()
    
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    scanner = CorrectedMosdepthRDScanner(args.coverage_file)
    scanner.coverage_threshold_factor = args.coverage_threshold
    scanner.known_rd_threshold = args.known_threshold
    scanner.min_rd_size = args.min_size
    scanner.max_rd_size = args.max_size
    scanner.merge_distance = args.merge_distance
    
    novel_rds, known_rds_analysis = scanner.run(args.known_rds)
    
    # Save results
    if not novel_rds.empty:
        novel_rds_formatted = novel_rds[[
            'chrom', 'start', 'end', 'size', 'type', 
            'coverage', 'chrom_coverage', 'proportion', 
            'deletion_category', 'known_rd'
        ]].copy()
        novel_rds_formatted.columns = [
            'Chromosome', 'Start', 'End', 'Size', 'Type', 
            'Coverage', 'ChromCoverage', 'Proportion', 
            'DeletionCategory', 'KnownRD'
        ]
        novel_rds_formatted.to_csv(args.output_novel, sep='\t', index=False)
    else:
        pd.DataFrame().to_csv(args.output_novel, sep='\t', index=False)
    
    if not known_rds_analysis.empty:
        known_rds_analysis.to_csv(args.output_known, sep='\t', index=False)
    else:
        pd.DataFrame().to_csv(args.output_known, sep='\t', index=False)
    
    # Print summary
    print(f"\nAnalysis completed")
    print(f"Found {len(novel_rds)} novel RD regions")
    
    if not known_rds_analysis.empty:
        absent_rds = known_rds_analysis[known_rds_analysis['Status'] == 'Absent']
        print(f"Known RDs: {len(absent_rds)} deletions confirmed, {len(known_rds_analysis) - len(absent_rds)} regions intact")
        print(f"\nConfirmed deletions: {', '.join(absent_rds['RD'].tolist())}")
    
    # Show validation statistics
    if not novel_rds.empty:
        print("\nNovel RD categories:")
        category_counts = novel_rds['deletion_category'].value_counts()
        for category, count in category_counts.items():
            print(f"  {category}: {count}")
        
        print("\nTop 5 strongest deletions:")
        top_deletions = novel_rds.nsmallest(5, 'proportion')[[
            'chrom', 'start', 'end', 'size', 'coverage', 'proportion', 'deletion_category', 'known_rd'
        ]]
        for _, row in top_deletions.iterrows():
            print(f"  {row['chrom']}:{row['start']}-{row['end']} - coverage: {row['coverage']:.2f}, proportion: {row['proportion']:.4f}")


if __name__ == '__main__':
    main()